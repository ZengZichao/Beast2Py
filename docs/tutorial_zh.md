# 教程

## 1. 安装

```bash
# 克隆仓库
git clone https://github.com/ZengZichao/Beast2Py.git
cd Beast2Py

# 安装
pip install -e .
```

未安装时，在仓库根目录下，下文每条命令都可以改用 `python3 -m beast2py.main <子命令>` 运行，入口完全相同。

## 2. 基本用法：单 partition 分析

### 步骤 1：准备序列数据

创建包含序列的 FASTA 文件：

```
>human
ATGGCAAATCATT...
>chimp
ATGGCAAATCATT...
>gorilla
ATGGCAAATCATT...
```

### 步骤 2：创建 YAML 配置

创建 `config.yaml`：

```yaml
metadata:
  analysis_name: "my_analysis"
  beast2_version: "2.7.8"

alignments:
  - id: "alignment"
    file: "input.fasta"
    data_type: "nucleotide"

partitions:
  - id: "alignment"
    site_model:
      substitution_model:
        type: "hky"
      gamma_categories: 4
      gamma_shape: {value: 0.5, lower: 0.0}
    clock_model:
      type: "strict"
    tree: "shared"

tree_prior:
  type: "yule"
  birth_rate: {value: 1.0, lower: 0.0}

calibrations:
  - name: "rootCalibration"
    taxa: null
    monophyletic: true
    distribution:
      type: "normal"
      parameters: {mean: 10.0, sigma: 1.0}
    provenance:
      source: "fossil"
      reference: "DOI:10.xxxx"
      calibration_type: "soft"

mcmc:
  chain_length: 10000000

output:
  file: "output.xml"
  embed_fingerprint: true
```

这里有两个默认值值得注意。

- 本例没有写 `mcmc.store_every`。它由链长按 `max(1, chain_length // 10000)` 推出（链长 1,000 万时为 `1000`），于是一次运行大约得到十个 `.xml.state` 检查点，而不是只有最末一个。
- 写成 `estimate: false` 的参数由生成器内联到引用它的元素里，因此不进入 `<state>`，不获得算子，也不出现在轨迹中。

### 步骤 3：生成 XML

```bash
beast2py generate --config config.yaml --output output.xml
```

`generate` 在写出文件之前依次运行三道验证闸门：Beast2Py 的结构与语义检查、校准冲突检测，以及（加上 `--beast2-validate` 时）BEAST2 解析器对每个对象调用的 `initAndValidate()`。所有运行过的闸门都通过后才写文件，任一失败都以非零码退出。不给 `--force` 时，`generate` 不会覆盖已存在的 `output.xml`。加上 `--diagnose`、`--fingerprint`、`--methods` 可以在同一次调用中产出报告与可重复性制品。

### 步骤 4：运行校准诊断

```bash
beast2py diagnose --config config.yaml --report report.html --sensitivity --visualize
```

### 步骤 5：运行 BEAST2

```bash
beast -overwrite output.xml
```

随机数种子是 BEAST2 的命令行参数，而不是 XML 的一部分（BEAST2 的 `MCMC` 没有 `seed` 输入，写出一个会让模型无法解析）：

```bash
beast -overwrite -seed 1 output.xml
```

每 `store_every` 个已记录样本，BEAST2 写出一个 `output.xml.state` 检查点，这正是中断的运行可以续算的原因。按上面的默认值，一条 1,000 万代的链大约会得到十个检查点。传 `mcmc.store_every: -1` 可关闭检查点，此时验证器会警告该链无法续算。

## 3. 快速模式

对于简单的单 partition 分析，使用 `quick` 命令：

```bash
beast2py quick \
    --alignment input.fasta \
    --tree-prior yule \
    --subst-model hky \
    --clock-model strict \
    --calibration-yaml calibrations.yaml \
    --chain-length 10000000 \
    --output quick_output.xml
```

这里的输出文件命名为 `quick_output.xml`，而不是沿用第 2 节的 `output.xml`：两个命令在没有 `--force` 时都不会
覆盖已存在的输出文件，因此复用同一个名字只会让这一步以拒绝写入收场，而不是把两次分析混进同一个文件。

此外，`quick` 还接受 `--gamma-categories N`（默认 0，表示不使用 Gamma）、`--pre-burnin N`、`--name NAME`、`--force` 与 `-v/--verbose`。`quick` 同样运行那三道闸门，`store_every` 的默认值也相同。

校准 YAML 文件：

```yaml
- name: "humanChimpMRCA"
  taxa: ["human", "chimp"]
  monophyletic: true
  distribution:
    type: "normal"
    parameters: {mean: 6.0, sigma: 0.5}
```

## 4. 多 Partition 分析

```yaml
alignments:
  - id: "gene1"
    file: "gene1.fasta"
    data_type: "nucleotide"
  - id: "gene2"
    file: "gene2.fasta"
    data_type: "nucleotide"

partitions:
  - id: "gene1"
    site_model:
      substitution_model: {type: "hky"}
    clock_model:
      type: "ucln"
      ucld_mean: {value: 1.0}
      ucld_stdev: {value: 0.333}
    tree: "shared"
  - id: "gene2"
    site_model:
      substitution_model: {type: "gtr"}
    clock_model:
      linked_to: "gene1"    # 共享 gene1 的时钟
    tree: "shared"
```

`tree: "shared"` 是唯一合法取值。逐分区的独立拓扑尚未实现，解析阶段会直接拒绝 `tree: separate`（以及在该位置写一个 partition id），而不是静默连锁。否则，本意为 *BEAST 风格的多基因座分析就变成了完全连锁的分析。在逐分区树实现之前，请为每个独立分区各运行一次分析。`linked_to` 必须指向另一个已存在的 partition，自链接与未知 id 都是错误。

## 5. 密码子分区

```yaml
alignments:
  - id: "alignment"
    file: "coding.fasta"
    data_type: "nucleotide"
  - id: "codon12"
    data_id: "alignment"
    filter: "1::3,2::3"
    data_type: "nucleotide"
  - id: "codon3"
    data_id: "alignment"
    filter: "3::3"
    data_type: "nucleotide"

partitions:
  - id: "codon12"
    site_model: {substitution_model: {type: "hky"}}
    clock_model: {type: "strict"}
    tree: "shared"
  - id: "codon3"
    site_model: {substitution_model: {type: "gtr"}}
    clock_model: {linked_to: "codon12"}
    tree: "shared"
```

声明 `filter` 时，`data_id` 为必填，且必须指向一个已声明的比对，不能指向过滤器自身。

Filter 表达式遵循 BEAST2 自己的 `FilteredAlignment.parseFilterSpec` 语法，坐标是 **1-based**（从 1 开始计数）。表达式是用逗号分隔的列表，每一项可以是：

- 单个位点（`N`）
- 闭区间（`from-to`，末端留空表示最后一个位点）
- 带步长的区间（`from-to\step`）
- 迭代器（`from:to:step`，上面两个例子用的就是 `from::step`）

于是 `1::3,2::3` 保留第 1、2、4、5、7、8 等位（每个密码子的前两位），`3::3` 保留 3、6、9 等位（第三位）。

解析器会针对源比对逐项校验范围。以下情况都会报错，并指名违规项：无法解析的表达式、步长小于 1、坐标超出比对长度、不选出任何位点。

## 6. 生成方法学描述

```bash
beast2py methods --config config.yaml --output methods.tex
```

此命令生成一段描述分析配置的 LaTeX 文本，适合在论文中使用。省略 `--output` 则改为直接打印。

## 7. 生成管道文件

```bash
# Snakemake 管道
beast2py pipeline --config config.yaml --type snakemake --output-dir pipeline

# Nextflow 管道
beast2py pipeline --config config.yaml --type nextflow --output-dir pipeline

# 两者都生成
beast2py pipeline --config config.yaml --type both --output-dir pipeline
```

`--output-dir` 默认为当前目录；`--type` 默认为 `snakemake`。

## 8. 验证

```bash
# 结构验证
beast2py validate --xml output.xml

# 同时用 BEAST2 验证
beast2py validate --xml output.xml --beast2
```

`--beast2` 运行无头的 BEAST2 解析器。它解析文件**并对每个对象调用 `initAndValidate()`**，因此检查的是模型能否初始化，而不只是能否解析。这与 `generate` 的闸门 3（`--beast2-validate`）走同一条代码路径，需要 BEAST2 以及 **JDK 17 或更新版本**。

随包发布的 `beast2_validate.sh`（位于包内，仓库根只有一个薄包装）按以下顺序查找 `BEAST.base.jar`：

1. `$BEAST2_JAR`
2. 最新的 `~/.beast/2.7/BEAST.base/*/lib/BEAST.base.jar`
3. 脚本旁边的 `tools/lib/` 或 `lib/`

脚本还会把 `~/.beast/2.7` 下安装的任何插件包追加到 class path。

脚本按一批候选路径搜索 JDK 17 或更新版本，而不是轻信 `PATH`——`PATH` 上一个古老的 Oracle `java` 垫片不该再遮蔽可用的 `openjdk@17`。搜索顺序为：`BEAST2_JAVA_CANDIDATES`（冒号分隔列表）、`JAVA_HOME`、`/usr/libexec/java_home`、Homebrew 与系统 JVM 目录，最后才是 `PATH` 上的 `java`。

两个开关可以移除对应的来源：`BEAST2_SKIP_JAVA_HOME_TOOL` 去掉 `/usr/libexec/java_home`，`BEAST2_SKIP_JDK_PROBES` 去掉 Homebrew 与系统 JVM 目录。设 `BEAST2_VALIDATE_TRACE` 会打印最终选中的 `java`。

退出码：`0` 表示全部有效，`1` 表示至少一个无效，`2` 表示环境搭建错误。脚本支持一次校验多个文件：

```bash
bash beast2_validate.sh output.xml examples/output_basic.xml
```

## 9. 列出支持的模型

```bash
beast2py list-models
```
