# Beast2Py 使用手册

- **版本：** 0.1.0
- **作者：** 曾子超 (Zichao Zeng) · [ORCID: 0000-0001-6553-970X](https://orcid.org/0000-0001-6553-970X)
- **日期：** 2026 年

---

## 目录

1. [引言](#1-引言)
2. [安装](#2-安装)
3. [两种调用模式：CLI 和 API](#3-两种调用模式cli-和-api)
4. [快速开始](#4-快速开始)
5. [CLI 命令行命令详解](#5-cli-命令行命令详解)
6. [Python API 详解](#6-python-api-详解)
7. [YAML 配置文件完整参考](#7-yaml-配置文件完整参考)
8. [替换模型](#8-替换模型)
9. [时钟模型](#9-时钟模型)
10. [树先验](#10-树先验)
11. [校准分布](#11-校准分布)
12. [多 Partition 分析](#12-多-partition-分析)
13. [末端日期](#13-末端日期)
14. [校准诊断](#14-校准诊断)
15. [可重复性工具](#15-可重复性工具)
16. [验证](#16-验证)
17. [MCMC 设置](#17-mcmc-设置)
18. [日志器配置](#18-日志器配置)
19. [初始化策略](#19-初始化策略)
20. [操作器权重](#20-操作器权重)
21. [示例配置文件](#21-示例配置文件)
22. [故障排除](#22-故障排除)
23. [附录：模型完整列表](#23-附录模型完整列表)

---

## 1. 引言

Beast2Py 是一个 Python 框架，用于生成、验证和诊断 BEAST2 分歧时间估计的 XML 配置。它提供全面的模型支持、自动化校准冲突检测、先验敏感性分析和可重复性工具。Beast2Py 同时提供命令行接口（CLI）和应用程序编程接口（API）两种调用模式，满足不同使用场景的需求。

### 核心功能

Beast2Py 从 YAML 配置生成完整的 BEAST2 XML，支持的模型包括：

- 13 种替换模型：7 种核苷酸模型和 6 种氨基酸模型
- 4 种时钟模型：严格时钟、UCLN、UCE、RLC
- 8 种树先验：Yule、CalibratedYule、Birth-Death、Coalescent 常量/指数、Bayesian Skyline、EBSP、BD Skyline Serial

该工具支持多 partition 分析（含密码子分区），并可在各 partition 之间链接或解链模型（link/unlink）。

诊断方面，Beast2Py 提供自动化校准冲突检测（时间一致性、分布重叠、单系性冲突）和先验敏感性分析（留一法 XML 生成）。此外，工具还支持校准分布可视化、分析指纹、方法学描述自动生成，以及 Snakemake/Nextflow 管道文件生成。

以下两条性质贯穿本手册的所有小节。

- **解析器失败即关闭。** 每个键都对照 schema 校验，每个布尔值都必须含义明确，每条比对都会检查等长、名称唯一和字符表合法，每个分布与树先验参数都会做值域检查。生成器无法兑现的输入一律报错，不会静默取默认值（见[第 22.1 节](#221-严格校验被拒绝的输入及其错误信息)）。
- **写出 XML 有闸门把关。** `generate` 与 `quick` 依次运行结构检查、校准冲突检测，以及可选的真实 BEAST2 解析与初始化检查。任一运行过的闸门失败，命令都以非零码退出，并且不写出文件（见[第 16 节](#16-验证)）。

### 设计理念

Beast2Py 采用分层架构设计，分离为 CLI/API 入口层、配置管理层、领域模块、数据模型层和工具支持层。领域模块中的诊断模块与可重复性模块是核心创新。模型注册表模式（`ModelRegistry`）把配置中的模型名称映射到对应的 BEAST2 XML spec 值，用户无需修改核心代码即可注册自定义模型。

### 与其他工具的关系

| 工具 | 语言 | 适用场景 |
|---|---|---|
| BEAUti2（BEAST2 自带） | 图形界面 | 交互式探索分析；v2.7+ 还内置交互式方法学段落查看器与先验分布预览图 |
| babette 套件 | R | R 一站式工作流：beautier 生成 XML、beastier 运行 BEAST2、tracerer 解析日志（ESS/HPD） |
| BEASTling / beast2-xml | Python | 配置文件/命令行式 XML 生成（面向语言学数据；FASTA/FASTQ 输入） |
| LPhy + LPhyBeast | Java | 概率模型规约语言 → BEAST XML（已面向 BEAST 3） |
| Beast2Py | Python | 声明式 YAML + 校准冲突检测 + 留一法先验敏感性 XML + 校准来源/指纹/方法学/管道制品 |

注意：仅先验采样并非 Beast2Py 独有——对任意既有 XML 可直接运行 `beast -sampleFromPrior output.xml`。Beast2Py 的贡献在于自动完成 N+1 留一法编排，并将冲突检测与可重复性制品整合到同一配置阶段工作流中。

---

## 2. 安装

### 2.1 前提条件

Beast2Py 需要 Python 3.10 或更高版本以及 pip 包管理器。生成、诊断与结构验证器均为纯 Python 实现。可选的 BEAST2 验证闸门还需要安装 BEAST2，以及 **JDK 17 或更新版本**的运行时，因为 BEAST 2.7.x 的 class 文件无法在更旧的 JDK 上加载。本文所有内容均已针对 **BEAST2 v2.7.8** 核实；其他 2.7.x 版本预期可用，但未经测试。

### 2.2 从源码安装

```bash
# 克隆仓库
git clone https://github.com/ZengZichao/Beast2Py.git
cd Beast2Py

# 以开发模式安装
pip install -e .
```

### 2.3 依赖包

以下依赖包将在安装时自动获取：

| 包 | 最低版本 | 用途 |
|---|---|---|
| pyyaml | >=6.0 | YAML 配置解析 |
| biopython | >=1.80 | 序列文件解析（FASTA/NEXUS） |
| numpy | >=1.20 | 数值计算和数组操作 |
| scipy | >=1.8 | 统计分布和 Wasserstein 距离计算 |
| matplotlib | >=3.5 | 校准分布可视化 |
| jinja2 | >=3.0 | HTML 诊断报告生成 |

### 2.4 可选依赖

| 包 | 用途 |
|---|---|
| lxml | 可选的 XML 解析后端，提供更快的解析速度 |
| pytest | 运行测试套件 |
| pytest-cov | 测试覆盖率统计 |
| black | 代码格式化 |
| flake8 | 代码风格检查 |
| mypy | 类型检查 |

安装开发依赖：

```bash
pip install -e ".[dev]"
```

安装可选 XML 后端：

```bash
pip install -e ".[xml]"
```

### 2.5 验证安装

```bash
# 检查版本
beast2py --version
# 输出: Beast2Py v0.1.0

# 列出支持的模型
beast2py list-models

# 验证 Python API
python3 -c "import beast2py; print(beast2py.__version__)"
# 输出: 0.1.0
```

未安装时，在克隆出的仓库根目录下，所有命令都可以改用 `python3 -m beast2py.main <子命令>` 运行。两种形式走的是同一个入口（`[project.scripts] beast2py = "beast2py.main:main"`）。

### 2.6 安装 BEAST2（可选，用于 BEAST2 验证闸门）

BEAST2 可从 https://www.beast2.org/ 下载。传入 `--beast2-validate`（用于 `generate` 与 `quick`）或 `--beast2`（用于 `validate`）即可启用第三道闸门；此时 `--beast2-path` 接受一个 BEAST2 的 `beast` 可执行文件，或一个 `beast2_validate.sh` 脚本。两者都不指定时，工具会自动定位随包发布的 `beast2_validate.sh`，依次查找安装包旁、项目根和 `PATH`。该脚本需要 `BEAST.base.jar`，按以下顺序查找：

1. `$BEAST2_JAR` 指向的路径（若设置了该变量）；
2. 最新的 `~/.beast/2.7/BEAST.base/*/lib/BEAST.base.jar`；
3. 脚本旁边的 `tools/lib/BEAST.base.jar` 或 `lib/BEAST.base.jar`。

脚本会把 `~/.beast/2.7` 下安装的插件包（BDSKY、NS 等）自动追加到 class path，因此依赖插件的那两个示例同样可以校验。

脚本需要 JDK 17 或更新版本，并会在一批**候选路径**中搜索，而不是轻信 `PATH`。原因是 macOS 机器常常在 `PATH` 上留有一个陈旧的 Oracle `java` 垫片，而可用的 `openjdk@17` 只是没有链接到 `PATH`。第一个主版本号不低于 17 的候选胜出，顺序为：

1. `BEAST2_JAVA_CANDIDATES` 中的条目（冒号分隔的列表，需显式设置）；
2. `$JAVA_HOME/bin/java`；
3. `/usr/libexec/java_home -v 17+` 报出的路径（可用 `BEAST2_SKIP_JAVA_HOME_TOOL=1` 跳过）；
4. Homebrew、`/usr/lib/jvm` 与 `/Library/Java/JavaVirtualMachines` 目录（可用 `BEAST2_SKIP_JDK_PROBES=1` 跳过）；
5. `PATH` 上的 `java`——排在最后，因为它是所有候选中最不明确的一个。

这些候选中没有一个达到 17 时，脚本才会放弃，以 **2** 退出，并列出未通过的候选：

```bash
bash beast2_validate.sh examples/output_basic.xml
# ... 以如下行结束：VALID: examples/output_basic.xml

BEAST2_VALIDATE_TRACE=1 bash beast2_validate.sh examples/output_basic.xml
# ... 在 stderr 上打印：[beast2_validate] using java 17 at /opt/homebrew/opt/openjdk@17/...
```

`JAVA_HOME` 仍然是强制指定 JDK 的手段；只有 JDK 8 的机器依旧会得到清晰的报错：

```bash
JAVA_HOME=/opt/homebrew/opt/openjdk@17/libexec/openjdk.jdk/Contents/Home \
    bash beast2_validate.sh examples/output_basic.xml
```

退出码为 `0`（全部有效）、`1`（至少一个无效）与 `2`（环境搭建错误：jar 或 JDK 缺失、版本过旧）。一次调用可以校验多个文件。

---

## 3. 两种调用模式：CLI 和 API

Beast2Py 提供两种等价的调用模式，用户可根据使用场景选择。

### 3.1 CLI 模式（命令行接口）

CLI 模式适用于命令行环境下的批量分析、脚本化工作流和管道集成。用户通过 `beast2py` 命令调用它，CLI 提供 8 个子命令：`generate`、`diagnose`、`quick`、`validate`、`methods`、`pipeline`、`fingerprint` 和 `list-models`。

```bash
beast2py generate --config config.yaml --output output.xml
```

### 3.2 API 模式（Python 编程接口）

API 模式适用于 Jupyter Notebook 交互式分析、Python 脚本中的程序化调用，以及与其他生物信息学工具的集成。导入 `beast2py` 包即可使用。

```python
from beast2py.api import Beast2Py

# 创建 API 实例
b2p = Beast2Py()

# 生成 XML
xml_str = b2p.generate_xml("config.yaml", output="output.xml")

# 运行诊断
report = b2p.diagnose("config.yaml", report_path="report.html")

# 验证 XML
result = b2p.validate_xml("output.xml")
```

也支持模块级便捷函数：

```python
from beast2py.api import generate_xml, diagnose, validate_xml

xml = generate_xml("config.yaml", output="output.xml")
report = diagnose("config.yaml", report_path="report.html")
result = validate_xml("output.xml")
```

### 3.3 两种模式的对应关系

| CLI 命令 | API 方法 | 模块级函数 |
|---|---|---|
| `beast2py generate` | `Beast2Py.generate_xml()` | `generate_xml()` |
| `beast2py quick` | `Beast2Py.quick_generate()` | `quick_generate()` |
| `beast2py diagnose` | `Beast2Py.diagnose()` | `diagnose()` |
| `beast2py validate` | `Beast2Py.validate_xml()` | `validate_xml()` |
| `beast2py methods` | `Beast2Py.generate_methods()` | `generate_methods()` |
| `beast2py pipeline` | `Beast2Py.generate_pipeline()` | `generate_pipeline()` |
| `beast2py fingerprint` | `Beast2Py.generate_fingerprint()` | `generate_fingerprint()` |
| `beast2py list-models` | `Beast2Py.list_models()` | `list_models()` |

---

## 4. 快速开始

### 4.1 基本工作流

Beast2Py 的典型工作流包括以下步骤。

1. 准备序列数据（FASTA 或 NEXUS 格式）。
2. 创建 YAML 配置文件，指定模型、校准和 MCMC 参数。
3. 用 Beast2Py 生成 BEAST2 XML。
4. 可选：运行校准诊断，检测潜在问题。
5. 用生成的 XML 运行 BEAST2，并分析结果。

### 4.2 CLI 快速开始

```bash
# 步骤 1：创建 YAML 配置文件（见第 7 节完整参考）
# 步骤 2：生成 XML（先过三道闸门；任一失败即不写出文件）
beast2py generate --config config.yaml --output output.xml -v

# 步骤 3：运行诊断（可选）
beast2py diagnose --config config.yaml --report report.html --sensitivity --visualize

# 步骤 4：验证 XML（可选，用 BEAST2 需 JDK 17）
beast2py validate --xml output.xml --beast2

# 步骤 5：运行 BEAST2
beast -overwrite output.xml
```

若 `output.xml` 已存在，`generate` 会直接拒绝并以退出码 1 结束，除非传入 `--force`。原因是同名的 `.fingerprint.json` / `.methods.tex` 仍然存在时，覆盖输出会让两次分析共用同一个名字。

在克隆本仓库后可直接运行、已实测的输出示例：

```text
$ beast2py generate -c examples/config_basic.yaml -o /tmp/output.xml -v
[INFO] Reading configuration: examples/config_basic.yaml
  [verbose] Analysis: primates_basic_calibration
  [verbose] Partitions: 1
  [verbose] Calibrations: 1
  [verbose] Tree prior: yule
[INFO] Generating BEAST2 XML...
[OK] Gate 1/3 passed: Beast2Py structural checks (ids, idrefs, priors, operators).
[OK] Gate 2/3 passed: no contradictory calibration conflicts.
[OK] XML written to: /tmp/output.xml
[OK] Done!
```

加 `--beast2-validate` 会多出一行 `[OK] Gate 3/3 passed: BEAST2 parsed and initialised the model.`。缺 BEAST2 时输出改为 `[ERROR] Gate 3/3 unavailable: …`，退出码为 **2**，并且不写任何文件。加上 `--allow-unvalidated` 会把该结果降级为 `[WARN]` 并照样写出，此时不应把该 XML 称为已通过 BEAST2 验证。

### 4.3 API 快速开始

```python
from beast2py.api import Beast2Py

b2p = Beast2Py()

# 一步完成生成、验证、指纹和方法学描述
result = b2p.generate(
    "config.yaml",
    output="output.xml",
    diagnose=True,
    fingerprint=True,
    methods=True,
)

print(f"XML: {result['xml_path']}")
print(f"验证通过: {result['validation'].is_valid}")
print(f"指纹: {result['fingerprint']}")
print(f"方法学: {result['methods_path']}")

# 查看诊断结果
if 'diagnostic_report' in result:
    report = result['diagnostic_report']
    for c in report.conflicts:
        print(f"[{c.severity}] {c.description}")
```

### 4.4 最小配置示例

```yaml
metadata:
  analysis_name: "example"

alignments:
  - id: "alignment"
    file: "primates.fasta"
    data_type: "nucleotide"

partitions:
  - id: "alignment"
    site_model:
      substitution_model: {type: "hky"}
    clock_model: {type: "strict"}
    tree: "shared"

tree_prior:
  type: "yule"
  birth_rate: {value: 1.0}

calibrations:
  - name: "root"
    taxa: null
    monophyletic: true
    distribution:
      type: "normal"
      parameters: {mean: 10.0, sigma: 1.0}

mcmc:
  chain_length: 10000000

output:
  file: "output.xml"
```

---

## 5. CLI 命令行命令详解

### 5.1 `generate` — 生成 BEAST2 XML

从 YAML 配置文件生成完整的 BEAST2 XML，可同时运行诊断、生成指纹和方法学描述。只有运行过的每道闸门全部通过，`generate` 才写出文件。

```bash
beast2py generate \
    --config config.yaml \
    --output output.xml \
    [--diagnose] \
    [--fingerprint] \
    [--methods] \
    [--beast2-validate] \
    [--beast2-path PATH] \
    [--force] \
    [--verbose]
```

**参数说明：**

| 参数 | 简写 | 必需 | 默认值 | 描述 |
|---|---|---|---|---|
| `--config` | `-c` | 是 | — | YAML 配置文件路径 |
| `--output` | `-o` | 是 | — | 输出 XML 文件路径 |
| `--diagnose` | | 否 | False | 同时运行校准诊断 |
| `--fingerprint` | | 否 | False | 生成分析指纹文件（`.fingerprint.json`） |
| `--methods` | | 否 | False | 生成 LaTeX 方法学描述（`.methods.tex`） |
| `--beast2-validate` | | 否 | False | 运行第三道闸门：BEAST2 解析器，**并对每个对象调用 `initAndValidate()`**（需 BEAST2 与 JDK 17） |
| `--beast2-path` | | 否 | `beast` | BEAST2 校验脚本（`.sh`）或可执行文件路径；缺省时自动定位随包发布的 `beast2_validate.sh`（项目根、安装包旁或 `PATH` 上） |
| `--force` | | 否 | False | `--output` 已存在时覆盖（默认：拒绝并退出码 1） |
| `--verbose` | `-v` | 否 | False | 显示详细输出信息 |

**闸门与退出码。** 闸门 1（Beast2Py 结构与语义自检）、闸门 2（校准冲突检测），以及按需开启的闸门 3（BEAST2 解析与初始化），全部在写文件**之前**运行。任一闸门失败都会打印 `[ERROR] Gate n/3 failed: …`，以退出码 1 返回，并且不写出文件；Python API 此时抛出 `ValueError`。若输出旁残留早先的 `.fingerprint.json` 或 `.methods.tex`，而本次调用省略了 `--fingerprint`/`--methods`，CLI 会警告这些文件正在变陈旧。三道闸门各自能证明什么、不能证明什么，见[第 16 节](#16-验证)。

**示例：**

```bash
# 基本生成
beast2py generate -c examples/config_basic.yaml -o /tmp/output.xml

# 生成并同时运行诊断、生成指纹和方法学描述
beast2py generate -c examples/config_basic.yaml -o /tmp/output.xml \
    --diagnose --fingerprint --methods -v

# 生成并用 BEAST2 原生验证（第三道闸门；需 JDK 17）
JAVA_HOME=/opt/homebrew/opt/openjdk@17/libexec/openjdk.jdk/Contents/Home \
    beast2py generate -c examples/config_basic.yaml -o /tmp/output.xml --beast2-validate

# 替换已存在的输出
beast2py generate -c examples/config_basic.yaml -o /tmp/output.xml --force
```

### 5.2 `diagnose` — 运行校准诊断

运行校准冲突检测、先验敏感性分析和分布可视化，生成自包含 HTML 报告。

```bash
beast2py diagnose \
    --config config.yaml \
    --report report.html \
    [--sensitivity] \
    [--visualize] \
    [--sensitivity-xml-dir DIR] \
    [--verbose]
```

**参数说明：**

| 参数 | 简写 | 必需 | 默认值 | 描述 |
|---|---|---|---|---|
| `--config` | `-c` | 是 | — | YAML 配置文件路径 |
| `--report` | | 是 | — | 输出 HTML 报告路径 |
| `--sensitivity` | | 否 | False | 生成先验敏感性分析 XML（N+1 个） |
| `--visualize` | | 否 | False | 生成校准分布密度图 |
| `--sensitivity-xml-dir` | | 否 | `<报告目录>/diagnostics` | 敏感性 XML 写入该目录下的 `sensitivity/` 子目录 |
| `--verbose` | `-v` | 否 | False | 显示详细输出信息 |

**示例：**

```bash
# 完整诊断
beast2py diagnose -c config.yaml --report report.html --sensitivity --visualize

# 仅冲突检测
beast2py diagnose -c config.yaml --report report.html
```

需要注意，`generate` 已经把冲突检测器作为它的闸门 2 运行过了；`diagnose` 提供的是独立而更完整的报告（敏感性 XML、密度图、HTML）。

### 5.3 `quick` — 快速模式

用户无需编写完整 YAML 配置，直接用命令行参数即可生成单 partition XML。`quick` 适用于简单分析和快速测试。

```bash
beast2py quick \
    --alignment input.fasta \
    --output output.xml \
    [--tree-prior yule] \
    [--subst-model hky] \
    [--clock-model strict] \
    [--gamma-categories 0] \
    [--calibration-yaml cals.yaml] \
    [--chain-length 10000000] \
    [--pre-burnin 0] \
    [--name "analysis_name"] \
    [--force] \
    [--verbose]
```

**参数说明：**

| 参数 | 简写 | 必需 | 默认值 | 描述 |
|---|---|---|---|---|
| `--alignment` | `-a` | 是 | — | 序列文件路径（FASTA/NEXUS），按核苷酸数据读取 |
| `--output` | `-o` | 是 | — | 输出 XML 文件路径 |
| `--tree-prior` | | 否 | `yule` | 树先验类型 |
| `--subst-model` | | 否 | `hky` | 替换模型类型 |
| `--clock-model` | | 否 | `strict` | 时钟模型类型 |
| `--gamma-categories` | | 否 | `0` | Gamma 速率类别数，0 表示不使用；大于 0 时自动创建初值为 0.5 的 `gamma_shape` 参数 |
| `--calibration-yaml` | | 否 | — | 校准点 YAML 文件路径 |
| `--chain-length` | | 否 | `10000000` | MCMC 链长度 |
| `--pre-burnin` | | 否 | `0` | MCMC 预烧入长度 |
| `--name` | | 否 | `quick_analysis` | 分析名称 |
| `--force` | | 否 | False | `--output` 已存在时覆盖 |
| `--verbose` | `-v` | 否 | False | 显示详细输出信息 |

`quick` 运行与 `generate` 相同的闸门 1 与闸门 2，同样遵守 `--force`，但它**没有 `--beast2-validate` 选项**。需要 BEAST2 校验时，改用 `beast2py validate --xml output.xml --beast2` 验证结果。

`quick` 按核苷酸数据读取比对，因此 `--subst-model` 只能指定核苷酸模型；数据类型交叉校验会拒绝氨基酸模型，这与 YAML 模式下的行为一致。同一个校验器也解析校准点 YAML，因此加引号的 `monophyletic: "false"`、反向的 `uniform` 界、未知分布类型在这里同样报错。

**校准点 YAML 文件格式（用于 `--calibration-yaml`）：**

```yaml
- name: "humanChimpMRCA"
  taxa: ["human", "chimp"]
  monophyletic: true
  distribution:
    type: "normal"
    parameters: {mean: 6.0, sigma: 0.5}
  provenance:
    source: "fossil"
    reference: "DOI:10.1038/nature12323"
    calibration_type: "soft"

- name: "rootCalibration"
  taxa: null
  monophyletic: true
  distribution:
    type: "uniform"
    parameters: {lower: 20.0, upper: 40.0}
  use_originate: true
```

**示例：**

```bash
# 最简单的快速生成
beast2py quick -a primates.fasta -o quick_output.xml

# 带校准和 Gamma
beast2py quick -a primates.fasta -o quick_full_output.xml \
    --tree-prior birth_death \
    --subst-model gtr \
    --clock-model ucln \
    --gamma-categories 4 \
    --calibration-yaml cals.yaml \
    --chain-length 50000000 \
    --name "my_analysis"
```

每个示例都写入各自的输出文件，因为 `quick` 与 `generate` 一样，在没有 `--force` 时不会替换已存在的输出文件。

### 5.4 `validate` — 验证 XML

先运行 Beast2Py 自身的结构与语义检查，再按需运行 BEAST2 解析与模型初始化检查。

```bash
beast2py validate \
    --xml output.xml \
    [--beast2] \
    [--beast2-path PATH]
```

**参数说明：**

| 参数 | 必需 | 默认值 | 描述 |
|---|---|---|---|
| `--xml` | 是 | — | 待验证的 XML 文件路径 |
| `--beast2` | 否 | False | 同时运行 BEAST2 解析 + 初始化检查（需 BEAST2 与 JDK 17） |
| `--beast2-path` | 否 | `beast` | BEAST2 可执行文件或 `beast2_validate.sh` 路径 |

**第 1 级（始终运行）。** 检查 XML 格式正确性与下列各项：

- 根元素必须为 `<beast>`
- 必需元素（`<data>` 与 `<run>`）存在
- id/idref/`@ref` 引用一致，无重复 ID
- `<state>` 元素存在，后验引用与日志器存在
- 序列数据非空且等长，分类单元名唯一

语义检查包括：每个被估计参数都带有 `Prior`，状态节点上下界未反向，初值落在自身界内且与声明维度匹配。`storeEvery` 关闭续算时，该级给出警告。

**第 2 级（`--beast2`）。** 无头 BEAST2 `XMLParser` 解析文件，**并对每个对象调用 `initAndValidate()`**。两级用不同措辞分别报告通过信息（`Level 1/2 passed: structural checks …` 与 `Level 2/2 passed: BEAST2 parsed and initialised the model.`）。任一级失败时，退出码为 1。

```bash
beast2py validate --xml examples/output_basic.xml
JAVA_HOME=/opt/homebrew/opt/openjdk@17/libexec/openjdk.jdk/Contents/Home \
    beast2py validate --xml examples/output_basic.xml --beast2
```

### 5.5 `methods` — 生成方法学描述

从配置自动生成适合论文使用的 LaTeX 方法学段落。

```bash
beast2py methods \
    --config config.yaml \
    [--output methods.tex]
```

| 参数 | 简写 | 必需 | 默认值 | 描述 |
|---|---|---|---|---|
| `--config` | `-c` | 是 | — | YAML 配置文件路径 |
| `--output` | `-o` | 否 | 标准输出 | 输出文件路径 |

**生成内容包括：** 替换模型描述（含引用）、时钟模型描述、树先验描述、校准点描述（含分布参数和来源）、MCMC 设置描述、分析指纹标识。

### 5.6 `pipeline` — 生成管道文件

生成 Snakemake 或 Nextflow 管道文件，实现端到端可重复分析工作流。

```bash
beast2py pipeline \
    --config config.yaml \
    --type snakemake|nextflow|both \
    [--output-dir .]
```

| 参数 | 必需 | 默认值 | 描述 |
|---|---|---|---|
| `--config` | 是 | — | YAML 配置文件路径 |
| `--type` | 否 | `snakemake` | 管道类型：`snakemake`、`nextflow` 或 `both` |
| `--output-dir` | 否 | `.` | 输出目录 |

管道包含以下规则或进程：XML 生成（调用 `beast2py generate`）、BEAST2 执行（调用 `beast -overwrite`）、校准诊断（调用 `beast2py diagnose`）、方法学描述生成（调用 `beast2py methods`）。

### 5.7 `fingerprint` — 生成分析指纹

`fingerprint` 生成基于配置 SHA-256 哈希的唯一分析标识符。

```bash
beast2py fingerprint \
    --config config.yaml \
    [--output fingerprint.json] \
    [--verbose]
```

| 参数 | 简写 | 必需 | 默认值 | 描述 |
|---|---|---|---|---|
| `--config` | `-c` | 是 | — | YAML 配置文件路径 |
| `--output` | `-o` | 否 | 仅输出到终端 | 输出 JSON 文件路径 |
| `--verbose` | `-v` | 否 | False | 显示详细指纹信息 |

**指纹格式：** `B2P-{配置哈希前 12 位}-{工具版本}`，例如 `B2P-ba18c26ed61f-0.1.0`。该标识符**不含日期与时区**，因此同一配置在同一比对数据上重跑，总能复现相同的字节与相同的标识符。12 个十六进制字符即 48 位摘要：按生日界估算，出现一次碰撞的概率达到 50% 需要 sqrt(2 * 2**48 * ln 2) = 1.98 × 10⁷ 个不同分析——即约 2 × 10⁷，而不是开方捷径给出的 2**24 = 1.7 × 10⁷。它是变更检测令牌，不是全局唯一键。

**指纹 JSON 文件内容**（对 `examples/config_basic.yaml` 的实测输出）：

```json
{
  "fingerprint": "B2P-04ab3d09277b-0.1.0",
  "config_hash": "04ab3d09277b4acb25c69b820e1f806b04adb48b08cd3af5ed9dc5d54a0eb5d1",
  "full_hash": "04ab3d09277b4acb25c69b820e1f806b04adb48b08cd3af5ed9dc5d54a0eb5d1",
  "data_hash": "30be5613d3ff6862",
  "hash_bits": 48,
  "collision_note": "The identifier truncates the configuration digest to 48 bits; a 50% chance of one collision needs sqrt(2 * 2**48 * ln 2), about 1.98e+07 distinct analyses, so it is a change-detection token rather than a globally unique key.",
  "tool_version": "0.1.0",
  "beast2_version": "2.7.8",
  "generation_time": "2026-09-29T18:01:17.890711",
  "analysis_name": "primates_basic_calibration",
  "xml_digest": "4f9352096f57259bcf987ba9853f943b406ee6c697f417185ed527c935ea0756"
}
```

上述 `generation_time` 是采集该样本的时刻；标识符本身不含时间成分，同一配置在同一比对数据上重跑总得到 `B2P-04ab3d09277b-0.1.0`。

`data_hash` 由每条比对内容的摘要再汇总而成，取前 16 位，会作为 `| Data: {hash16}` 附在 XML 的指纹注释里。`generation_time` **只**出现在这个侧车文件中；XML 本体带的是无时间戳的注释（见第 15.1 节）。

### 5.8 `list-models` — 列出支持的模型

显示 Beast2Py 支持的全部模型类型。

```bash
beast2py list-models
```

此命令无需参数，输出包括替换模型、时钟模型、树先验、校准分布和数据类型的完整列表。

### 5.9 全局选项

```bash
beast2py --version       # 显示版本号
beast2py --help          # 显示帮助信息
beast2py <command> --help  # 显示特定命令帮助
```

---

## 6. Python API 详解

### 6.1 Beast2Py 类

`Beast2Py` 类是 API 的核心入口，封装了所有功能。

```python
from beast2py.api import Beast2Py

b2p = Beast2Py(base_dir=".")  # base_dir 用于解析相对文件路径
```

**构造函数参数：**

| 参数 | 类型 | 默认值 | 描述 |
|---|---|---|---|
| `base_dir` | `str` | 当前工作目录 | 解析配置中相对文件路径的基准目录 |

### 6.2 `generate_xml` — 生成 XML

```python
xml_str = b2p.generate_xml(
    config,           # 配置文件路径、字典或 BEASTConfig 对象
    output=None,      # 可选：输出文件路径
    force=False,      # 可选：覆盖已存在的输出文件
)
```

`config` 参数接受三种类型：YAML 文件路径（字符串）、配置字典（`dict`）或已解析的 `BEASTConfig` 对象。返回生成的 XML 字符串。

给出 `output` 时，`generate_xml()` 让 XML 走**与 CLI 相同的闸门**（结构检查，随后是校准冲突检测），全部通过才写盘。任一闸门失败都会抛出 `ValueError`，并保持该文件原样；目标文件已存在时，需要 `force=True` 才会覆盖。这是刻意的设计：早期版本先写文件再检查，于是库与命令行对“已验证”的含义并不一致。如果只需要未经闸门的原始字符串，请直接调用 `XMLWriter(config).generate_xml()`。

### 6.3 `quick_generate` — 快速生成

```python
xml_str = b2p.quick_generate(
    alignment="primates.fasta",
    output="output.xml",
    tree_prior="yule",
    subst_model="hky",
    clock_model="strict",
    gamma_categories=0,
    calibration_yaml="cals.yaml",  # 可选
    chain_length=10_000_000,
    pre_burnin=0,
    name="my_analysis",
    force=False,
)
```

它运行同样的三道闸门，`force` 的语义也与 `generate_xml()` 一致。

### 6.4 `diagnose` — 运行诊断

```python
report = b2p.diagnose(
    config,
    output_dir="diagnostics",
    run_sensitivity=True,
    run_visualization=True,
    report_path="report.html",  # 可选
)

# report 是 DiagnosticReport 对象
for conflict in report.conflicts:
    print(f"[{conflict.severity}] {conflict.conflict_type}: {conflict.description}")

for summary in report.summaries:
    print(f"{summary.name}: mean={summary.stats['mean']:.2f}")

print(f"敏感性 XML: {report.sensitivity_xmls}")
print(f"可视化文件: {report.visualization_file}")
```

### 6.5 `detect_conflicts` — 仅检测冲突

```python
conflicts = b2p.detect_conflicts(config)
for c in conflicts:
    print(f"[{c.severity}] {c.calibration_a} ↔ {c.calibration_b}: {c.description}")
```

### 6.6 `generate_sensitivity_xmls` — 生成敏感性 XML

```python
xml_paths = b2p.generate_sensitivity_xmls(config, output_dir="sensitivity")
print(f"生成了 {len(xml_paths)} 个 XML 文件")
```

### 6.7 `plot_calibrations` — 绘制校准分布

```python
b2p.plot_calibrations(config, output_file="calibrations.png")
```

支持 PNG、PDF 和 SVG 格式。

### 6.8 `validate_xml` — 验证 XML

```python
# 验证文件
result = b2p.validate_xml("output.xml")

# 验证 XML 字符串
result = b2p.validate_xml(xml_string)

# 同时使用 BEAST2 原生验证
result = b2p.validate_xml("output.xml", beast2_validate=True, beast2_path="/opt/beast/bin/beast")

print(f"有效: {result.is_valid}")
print(f"错误: {result.errors}")
print(f"警告: {result.warnings}")
```

### 6.9 `generate_fingerprint` — 生成指纹

```python
fingerprint = b2p.generate_fingerprint(config, output="fingerprint.json")
print(f"指纹: {fingerprint}")
```

### 6.10 `generate_methods` — 生成方法学描述

```python
methods_text = b2p.generate_methods(config, output="methods.tex")
print(methods_text[:200])
```

### 6.11 `generate_pipeline` — 生成管道文件

```python
# 生成 Snakemake
path = b2p.generate_pipeline(config, pipeline_type="snakemake", output_dir="pipeline/")

# 生成 Nextflow
path = b2p.generate_pipeline(config, pipeline_type="nextflow", output_dir="pipeline/")

# 同时生成两者
sm_path, nf_path = b2p.generate_pipeline(config, pipeline_type="both", output_dir="pipeline/")
```

### 6.12 `generate` — 一站式生成

`generate` 方法是 `generate_xml` 的增强版，可同时完成验证、诊断、指纹和方法学描述。

```python
result = b2p.generate(
    config="config.yaml",
    output="output.xml",
    diagnose=True,
    fingerprint=True,
    methods=True,
    beast2_validate=False,
    beast2_path="beast",
)

# 返回字典包含：
# 'xml_path': XML 文件路径
# 'validation': ValidationResult 对象
# 'fingerprint': 指纹字符串（如果 fingerprint=True）
# 'fingerprint_path': 指纹文件路径（如果 fingerprint=True）
# 'methods_path': 方法学文件路径（如果 methods=True）
# 'diagnostic_report': DiagnosticReport 对象（如果 diagnose=True）
# 'beast2_validation': {'success': bool, 'output': str}（如果 beast2_validate=True）
```

`generate()` 把写盘动作委托给 `generate_xml()`，因此继承了同样的闸门与 `force` 行为：任一闸门失败都会先抛出 `ValueError`，不会写出文件。设 `beast2_validate=True` 而环境中没有 BEAST2 时，第三道闸门报告的是*跳过*而非通过；`XMLValidator.beast2_validation_status()` 区分 `"passed"`、`"failed"` 与 `"unavailable"`，因此缺失校验器时，该状态报告为不可用，绝不会报告为成功。

### 6.13 `list_models` — 列出模型

```python
models = b2p.list_models()
print(models["substitution_models"])  # ['blosum62', 'cprev', 'dayhoff', 'gtr', 'hky', ...]
print(models["clock_models"])          # ['rlc', 'strict', 'uce', 'ucln']
print(models["tree_priors"])           # ['bayesian_skyline', 'bd_skyline_serial', ...]
```

### 6.14 `parse_config` — 解析配置

```python
config = b2p.parse_config("config.yaml")
print(f"分析名称: {config.metadata['analysis_name']}")
print(f"Partition 数: {len(config.partitions)}")
print(f"校准点数: {len(config.calibrations)}")
print(f"树先验: {config.tree_prior_type.value}")
```

### 6.15 `read_sequence` — 读取序列

```python
alignment = b2p.read_sequence("primates.fasta", data_type="nucleotide")
print(f"分类单元数: {alignment.n_taxa}")
print(f"位点数: {alignment.n_sites}")
print(f"分类单元名: {alignment.taxa_names}")
```

### 6.16 模块级便捷函数

所有 `Beast2Py` 类方法都有对应的模块级函数，无需创建实例即可调用：

```python
from beast2py.api import (
    generate_xml,
    quick_generate,
    diagnose,
    validate_xml,
    generate_fingerprint,
    generate_methods,
    generate_pipeline,
    list_models,
    parse_config,
    read_sequence,
)

# 直接调用
xml = generate_xml("config.yaml", output="output.xml")
result = validate_xml("output.xml")
models = list_models()
```

### 6.17 直接使用配置类型

API 同时提供所有配置数据类，可用于程序化构建配置：

```python
from beast2py import (
    BEASTConfig,
    Partition,
    SiteModelConfig,
    ClockModelConfig,
    ClockModelType,
    TreePriorType,
    MCMCConfig,
    CalibrationPoint,
    DistributionConfig,
    RealParameter,
)

# 程序化构建配置
site_model = SiteModelConfig(
    substitution_model={"type": "hky"},
    gamma_categories=4,
    gamma_shape=RealParameter(value=0.5, lower=0.0),
)

clock_model = ClockModelConfig(type=ClockModelType.UCLN)

# ... 构建 Partition, Alignment 等对象
```

---

## 7. YAML 配置文件完整参考

### 7.1 完整配置模板

```yaml
# === 元数据 ===
metadata:
  analysis_name: "my_analysis"          # 分析名称（必填）
  analysis_description: "分析描述"       # 分析描述（可选）
  author: "你的名字"                     # 作者（可选）
  date: "2026-08-01"                    # 日期（可选）
  beast2_version: "2.7.8"              # BEAST2 版本（默认 2.7.8）
  tool_version: "0.1.0"                # Beast2Py 版本（默认 0.1.0）

# === 序列数据 ===
alignments:
  - id: "gene1"                         # 比对 ID（必填）
    file: "gene1.fasta"                 # 序列文件路径（必填，除非使用 filter）
    format: "auto"                      # 格式：fasta / nexus / auto（默认 auto）
    data_type: "nucleotide"             # 数据类型：nucleotide / aminoacid（默认 nucleotide）
    # 以下用于 FilteredAlignment（密码子分区）：
    filter: "3::3"                      # 过滤表达式（可选）
    data_id: "alignment"                # 基础比对 ID（使用 filter 时必填）

# === Partition 配置 ===
partitions:
  - id: "gene1"                         # Partition ID（必填，需与比对 ID 匹配）
    site_model:
      substitution_model:
        type: "hky"                     # 模型类型（见第 8 节）
        # 模型特定参数（见第 8 节）
      gamma_categories: 4               # Gamma 速率类别数（0=无，默认 0）
      gamma_shape:                      # Gamma 形状参数（gamma_categories>0 时）
        value: 0.5
        lower: 0.0
        estimate: true                  # 是否估计（默认 true）
      proportion_invariant:             # 不变位点比例（可选）
        value: 0.01
        lower: 0.0
        upper: 1.0
        estimate: true
    clock_model:
      type: "strict"                    # 时钟类型（见第 9 节）
      # 时钟特定参数（见第 9 节）
      linked_to: null                   # 链接到另一个 partition 的时钟（可选）
    tree: "shared"                      # "shared"=共享树（默认）

# === 末端日期（可选） ===
tip_dates:
  enabled: false                        # 是否启用（默认 false）
  trait_name: "date-forward"            # date-forward / date-backward
  units: "year"                         # year / month / day
  dates:                                # 分类单元到日期的映射
    taxonA: 2000
    taxonB: 2005
  # file: "tip_dates.csv"              # 替代：从 CSV 文件加载

# === 树先验 ===
tree_prior:
  type: "yule"                          # 树先验类型（见第 10 节）
  # 树先验特定参数（见第 10 节）

# === 校准点 ===
calibrations:
  - name: "humanChimpMRCA"              # 校准点名称（必填）
    taxa: ["human", "chimp"]            # 分类单元列表（null=根节点）
    monophyletic: true                  # 是否强制单系（默认 true）
    distribution:
      type: "normal"                    # 分布类型（见第 11 节）
      parameters: {mean: 6.0, sigma: 0.5}  # 分布参数
      offset: 0.0                       # 偏移量（默认 0.0）
      hyperpriors:                      # 超参数先验（可选）
        mean:
          type: "uniform"
          parameters: {lower: 0.0, upper: 20.0}
    use_originate: false                # 是否使用 stem 节点（默认 false）
    provenance:                         # 来源元数据（可选）
      source: "fossil"                  # 来源类型：fossil / biogeographic / secondary / other
      reference: "DOI:10.1038/xxxxx"    # 参考文献 DOI 或 URL
      calibration_type: "soft"          # soft / hard
      original_age_min: 5.5            # 原始最小年龄
      original_age_max: 7.0            # 原始最大年龄
      notes: "化石描述"

# === 校准方法 ===
calibration_method: "mrca_prior"        # mrca_prior / calibrated_yule（默认 mrca_prior）

# === 参数先验（可选） ===
parameter_priors:
  - parameter: "birthRate"              # 参数名
    distribution:
      type: "uniform"
      parameters: {lower: 0.0, upper: 1000.0}

# === MCMC 设置 ===
mcmc:
  type: "standard"                      # standard / nested_sampling（默认 standard）
  chain_length: 10000000                # 链长度（默认 10000000）
  pre_burnin: 0                         # 预烧入（默认 0，须 < chain_length）
  store_every: 10000                    # 检查点间隔，按"已记录样本数"计（缺省随链长伸缩 max(1, chain_length//10000)；-1 关闭）
  sample_from_prior: false              # 是否仅从先验采样（默认 false）
  seed: 1                               # 供 `beast -seed` 使用；**不写入 XML**（BEAST2 的 MCMC 无 seed 输入）
  # Nested Sampling 参数（type=nested_sampling 时）：
  particle_count: 1                     # 粒子数（默认 1）
  sub_chain_length: 10000               # 子链长度（默认 10000）

# === 日志器 ===
loggers:
  trace_log:
    file_name: "output.$(seed).log"     # 日志文件名（支持 $(seed) 变量，由 BEAST2 在运行时展开）
    log_every: 1000                     # 记录间隔
  tree_log:
    file_name: "output.$(seed).trees"
    log_every: 1000
  screen_log:
    log_every: 10000                    # 屏幕日志间隔

# === 初始化 ===
initialization:
  tree_type: "random"                   # random / upgma / newick（默认 random）
  newick_file: null                     # newick 类型时指定树文件

# === 操作器权重（可选覆盖） ===
operator_weights:
  tree_scaler: 3.0
  uniform: 30.0
  # ... 其他操作器

# === 校准诊断设置（可选，供闸门 2 使用） ===
diagnostics:
  overlap_measure: "relative"           # relative（默认，无量纲）/ absolute / both
  overlap_threshold: 2.0                # 绝对判据：Wasserstein 距离（时间单位）
  relative_overlap_threshold: 0.15      # 无量纲判据：W / 两先验 95% 区间的平均宽度

# === 输出设置 ===
output:
  file: "output.xml"                    # 输出 XML 文件路径
  embed_fingerprint: true               # 是否在 XML 中嵌入指纹注释（默认 true）
```

### 7.2 配置参数说明

**RealParameter 格式：** 许多参数接受 `RealParameter` 格式，可以是简单标量或完整字典：

```yaml
# 简单标量
clock_rate: 1.0

# 完整字典
clock_rate:
  value: 1.0        # 参数值
  lower: 0.0        # 下界（可选）
  upper: 10.0       # 上界（可选）
  dimension: 1      # 维度（默认 1）
  estimate: true    # 是否估计（默认 true）
```

**边界是覆盖项，而非必填项。** 不写 `lower`/`upper` 时，生成器按参数的统计角色补默认支撑集（`model.py::PARAM_SUPPORT`）：频率与概率 `[0, 1]`，速率与群体大小 `> 0`，计数与组大小 `>= 1`，`growth_rate` 等有符号量不设界。写出的 `lower`/`upper` 会替换默认值。解析阶段会拒绝两类输入：`value` 落在结果区间之外，`lower >= upper`。

**`estimate: false` 对每一类参数都真正生效** —— 替换模型速率、碱基频率、gamma 形状、不变位点比例、时钟速率、`ucld.mean` 与 `ucld.stdev`。生成器把固定参数以内联形式写出（`<parameter estimate="false">`），因此它不会成为状态节点，不会获得算子，也不会出现在 trace 里。如果又为该参数写了 `parameter_priors` 条目，该先验就指向一个不存在的目标，解析器会拒绝这次运行——这正是设计目的：把先验删掉。

```yaml
substitution_model:
  type: "hky"
  kappa: {value: 2.0, estimate: false}   # 固定：无 kappaScaler、无 trace 列
```

**布尔值** 必须是真正的 YAML 布尔，或 `true/false/yes/no/on/off/1/0` 之一。解析器接受加引号的形式，并按字面含义解析。Python 中 `bool("false")` 为 `True`，因此过去加引号的 `"false"` 会反转设置。现在 `"false"` 按 false 解析，解析器会拒绝 `"nope"` 这类含糊值。

---

## 8. 替换模型

### 8.1 核苷酸替换模型

表中的“参数”一列，列出每个模型在 `substitution_model:` 内**唯一**接受的键（除 `type` 之外，使用频率的模型还可带 `frequencies`）。任何其他键——例如在不容纳该键的模型上写 `rates`——解析器一律拒绝。

| 模型 | 配置名 | 参数 | 描述 |
|---|---|---|---|
| Jukes-Cantor | `jc69` | （无） | 等速率、等频率的最简模型；spec `JukesCantor` |
| HKY | `hky` | `kappa` | 转换/颠换偏倚模型 |
| TN93 | `tn93` | `kappa1`、`kappa2` | 两个转换速率 |
| GTR | `gtr` | `rateAC`、`rateAG`、`rateAT`、`rateCG`、`rateCT`、`rateGT`（或 `rates`） | 通用时间可逆模型 |
| SYM | `sym` | 六个 `rate*` 键（或 `rates`） | 对称模型：六个自由交换率、**碱基频率相等**——需显式指定，见下方注记 |
| TIM | `tim` | `rateAG`、`rateCT`、`rateTransversions1`、`rateTransversions2` | 转换模型（两个相等的颠换率 ×2） |
| TVM | `tvm` | `rateAC`、`rateAT`、`rateCG`、`rateGT`、`rateTransitions` | 颠换模型（四个自由颠换率 + 一个共享转换率） |

**`rates:` 仅是 GTR/SYM 的便捷写法。** `rates:` 必须按 `rateAC rateAG rateAT rateCG rateCT rateGT` 的顺序列出**恰好六个**值，且不能与显式命名的 `rate*` 键混用，否则展开时会静默覆盖其中一项。解析器不再为短向量补 1.0，也不再截断长向量。TIM 与 TVM 施加的是速率**等式约束**，六个自由速率向量无法表达，因此在这两个模型上写 `rates` 时，解析器会报错，并给出应改用的参数名。

```yaml
# GTR via the vector
substitution_model:
  type: "gtr"
  rates: {value: "1.0 4.0 1.0 1.0 4.0 1.0", lower: 0.0}
  frequencies: {mode: "estimated", value: 0.25, dimension: 4}

# SYM via the vector
substitution_model:
  type: "sym"
  rates: {value: "1.0 4.0 1.0 1.0 4.0 1.0", dimension: 6, lower: 0.0}

# TIM with its own parameters (see examples/config_subst_tim.yaml)
substitution_model:
  type: "tim"
  rateAG: {value: 4.0, lower: 0.0}
  rateCT: {value: 4.0, lower: 0.0}
  rateTransversions1: {value: 1.0, lower: 0.0}
  rateTransversions2: {value: 1.0, lower: 0.0}
  frequencies: {mode: "empirical"}

# TVM with its own parameters (see examples/config_subst_models.yaml)
substitution_model:
  type: "tvm"
  rateAC: {value: 1.0, lower: 0.0}
  rateAT: {value: 1.0, lower: 0.0}
  rateCG: {value: 1.0, lower: 0.0}
  rateGT: {value: 1.0, lower: 0.0}
  rateTransitions: {value: 4.0, lower: 0.0}
```

**平衡频率（`frequencies:`）。** 频率只适用于需要频率的核苷酸模型（`hky`、`tn93`、`gtr`、`sym`、`tim`、`tvm`）；`jc69` 与氨基酸矩阵自带频率。共有四种显式模式：

| `mode` | 产出的 XML | 含义 |
|---|---|---|
| `estimated`（默认） | `<frequencies frequencies="@…freqParameter">` + 状态节点 | 自由概率向量，归一化到和为 1，参与采样并记录 |
| `fixed` | `<frequencies><parameter estimate="false" value="…"/></frequencies>` | 生成器把用户向量原样写出，置于状态空间之外，既不采样也不记录 |
| `uniform` | `<frequencies spec="Frequencies" data="@aln" estimate="false"/>` | BEAST2 的“字符均匀分布”频率 |
| `empirical` | `<frequencies spec="Frequencies" data="@aln" estimate="true"/>` | BEAST2 依据比对计数得到的频率 |

> **为何引入 `mode`。** 在 BEAST2 中，`<frequencies estimate="false">` 的含义是*字符均匀分布*，**不**表示“用我给的数值”。因此过去 `estimate: false` 会静默地丢弃用户提供的向量。现在请明确表达意图：`mode: fixed` 写出你的向量，`mode: uniform` 向 BEAST2 索取等频率，`mode: empirical` 让 BEAST2 计数。为向后兼容，未显式给 `mode` 而在向量上写 `estimate: false` 时，解析器会把它读作 `mode: fixed`。`uniform` 与 `empirical` 需要 `data_id` 指向比对。`fixed`/`estimated` 向量必须是合法概率向量：逐项大于 0、长度匹配 `dimension`、和为 1（容差 2%，随后归一化）。

> **SYM 的注意事项。** SYM *就是*“碱基频率相等的 GTR”，但生成器不会主动替你施加这一约束：没有任何 `frequencies:` 块时，它产出的是一个自由的（被估计的）频率向量，于是模型实际上成了 GTR。若要得到 `examples/config_subst_sym.yaml` 那样的等频率 SYM，请写 `frequencies: {mode: uniform}`。

### 8.2 氨基酸替换模型

| 模型 | 配置名 | 描述 |
|---|---|---|
| WAG | `wag` | WAG 矩阵 |
| JTT | `jtt` | JTT 矩阵 |
| Dayhoff | `dayhoff` | Dayhoff 矩阵 |
| BLOSUM62 | `blosum62` | BLOSUM62 矩阵 |
| CPREV | `cprev` | CPREV 矩阵 |
| MTREV | `mtrev` | MTREV 矩阵（线粒体） |

```yaml
substitution_model: {type: "wag"}
```

这些都是自带经验平衡频率的 `EmpiricalSubstitutionModel` 实例，因此既不生成 `<frequencies>` 元素，也不接受任何 `frequencies:` 块。每个模型在解析阶段都做三方交叉校验：声明的 `data_type`、模型所属家族，以及比对中*实际观察到*的字符表。因此在核苷酸数据上写 `wag`、在氨基酸数据上写 `hky`、或对一份只含 `A C G T` 的文件声明 `data_type: aminoacid`，解析器都会报错。错误信息直指问题本身，而不是在离病因很远的地方抛出一个 BEAST2 解析错误。

```yaml
# Amino-acid data: the alignment must say so, and only amino-acid models are allowed
alignments:
  - id: "protein"
    file: "examples/aminoacid.fasta"     # 10 taxa x 234 aa
    data_type: "aminoacid"
partitions:
  - id: "protein"
    site_model:
      substitution_model: {type: "jtt"}
      gamma_categories: 4
      proportion_invariant: {value: 0.05, lower: 0.0, upper: 1.0}
```

### 8.3 位点模型选项

位点模型在替换模型之上添加速率异质性：

```yaml
site_model:
  substitution_model: {type: "hky"}
  gamma_categories: 4          # Gamma 速率类别数（0 = 无 Gamma）
  gamma_shape:                 # Gamma 形状参数
    value: 0.5
    lower: 0.0
    estimate: true             # 估计或固定（默认 true）
  proportion_invariant:        # 不变位点比例（可选）
    value: 0.01
    lower: 0.0
    upper: 1.0
    estimate: true
```

**Gamma+I 组合：** 同时设置 `gamma_categories > 0` 和 `proportion_invariant` 即可实现 Gamma+I 模型。

**`gamma_shape` 需要两个以上类别。** BEAST2 的 `SiteModel` 把形状输入注明为 `Ignored if gammaCategoryCount 1 or less`，因此在 `gamma_categories <= 1` 时给出 `gamma_shape`，解析器会拒绝。此时被估计的形状参数只是一条纯先验随机游走，对似然毫无贡献，却仍然在轨迹里占一列，容易被误读成“γ 形状参数”。

---

## 9. 时钟模型

| 模型 | 配置名 | 关键参数 | 描述 |
|---|---|---|---|
| 严格时钟 | `strict` | `clock_rate` | 所有分支相同速率 |
| UC 松弛对数正态 | `ucln` | `ucld_mean`、`ucld_stdev` | 不相关松弛对数正态时钟 |
| UC 松弛指数 | `uce` | `ucld_mean` | 不相关松弛指数时钟 |
| 随机局部时钟 | `rlc` | （自动） | 随机局部时钟模型 |

`ucld_mean` 与 `ucld_stdev` **仅**对 `ucln` 与 `uce` 有效；在 `strict` 或 `rlc` 时钟上写其中之一时，解析器会报错（`… only applies to the ucln/uce clocks; the strict clock has no such parameter`）。UCLN 的速率分布是 `M = 1` 固定、`S = ucld.stdev` 的对数正态；UCE 的速率分布是均值为 1 的指数分布。

**配置示例：**

```yaml
# 严格时钟
clock_model:
  type: "strict"
  clock_rate: {value: 1.0, lower: 0.0}

# UCLN 松弛时钟 —— 上下界都被兑现，ucld.stdev 因此无从漂移
clock_model:
  type: "ucln"
  ucld_mean: {value: 1.0, lower: 0.0, upper: 10.0}
  ucld_stdev: {value: 0.333, lower: 0.0, upper: 1.0}

# UCE 松弛指数时钟
clock_model:
  type: "uce"
  ucld_mean: {value: 1.0, lower: 0.0}

# 随机局部时钟
clock_model:
  type: "rlc"

# 链接到另一个 partition 的时钟
clock_model:
  linked_to: "gene1"    # 共享 gene1 的时钟模型
```

**注意**

- **生成器会兑现上下界。** 过去生成器会在 `clock_rate`、`ucld.mean` 与 `gamma_shape` 上丢弃 `upper`，使 `ucld.stdev` 上方无界——这正是经典的 UCLN 标准差漂移。如今 `lower`/`upper` 都会传入 XML；不给任何界时，生成器按参数角色套用默认界（`ucld.mean`、`clock.rate` 大于 0；`ucld.stdev` 落在 `[0, 1]`）。
- **固定的时钟速率不能是唯一的时间尺度。** `clock_rate: {estimate: false}` 如今同样得到兑现：该参数以 `estimate="false"` 内联写出，离开状态空间、算子集合与轨迹。如果此时没有给出任何校准，绝对时间不可识别，解析器会拒绝这次运行（错误信息见[第 22.1 节](#221-严格校验被拒绝的输入及其错误信息)）。
- **`ucld.mean` 归一化。** 松弛时钟 XML 以 `normalize="true"` 产出，与 BEAUti 的参考版式一致，从而保持速率与节点年龄可识别。

---

## 10. 树先验

| 树先验 | 配置名 | 关键参数 | 描述 |
|---|---|---|---|
| Yule | `yule` | `birth_rate` | 纯生过程 |
| CalibratedYule | `calibrated_yule` | `birth_rate` | 集成校准的 Yule 模型 |
| Birth-Death | `birth_death` | `birth_rate`、`death_rate`、`sampling_rate` | 生灭过程 |
| Coalescent 常量 | `coalescent_constant` | `pop_size` | 常量群体溯祖 |
| Coalescent 指数 | `coalescent_exponential` | `pop_size`、`growth_rate` | 指数增长溯祖 |
| Bayesian Skyline | `bayesian_skyline` | （自动计算） | 贝叶斯天际线图 |
| EBSP | `ebsp` | （自动计算） | 扩展贝叶斯天际线图（BEAST 2.7 核心包内置） |
| BD Skyline Serial | `bd_skyline_serial` | `birth_rate`、`death_rate`、`sampling_rate`、`rho` | 分段常量 Birth-Death 天际线，支持序列化取样（需 BDSKY 插件包） |

**允许域（admissible domains）。** 解析阶段会校验以下约束。若不校验，XML 照样写出，只在其后的 BEAST2 内部报错：

- `birth_death` 要求 `0 <= death_rate < birth_rate` **且** `0 < sample_probability <= 1`。`BirthDeathGernhard08Model` 以 `birth - death`（>= 0）与 `death / birth`（< 1）参数化，因此超临界过程（`death_rate >= birth_rate`）在该模型中没有合法状态。错误信息会指向 `bd_skyline_serial`（BDSKY），后者分别估计成种、灭绝与取样率，允许 `death >= birth`。`birth_rate` 与 `death_rate` 也必须一起估计或一起固定，因为只有派生量才作为状态节点存在。
- `bd_skyline_serial` 要求 `birth_rate >= 0`、`death_rate >= 0`，且 `sampling_rate` 与 `rho` 落在 `[0, 1]`。
- `yule`、`calibrated_yule`、`birth_death` 需要 `birth_rate > 0`；溯祖先验需要 `pop_size > 0`。Yule 类先验的 `birth_rate` 是先验自身的参数，对其写 `estimate: false` 时，解析器会报错，而不是静默忽略。
- `coalescent_exponential` 的 `growth_rate` 是*有符号*的（负值表示群体下降），因此不设界，并配一个宽正态默认先验。

**默认先验。** 若某个被估计的连续参数缺少 `parameter_priors` 条目，生成器会为其指定一个弱信息默认先验（`xml_writer.py::_default_prior_element`）：正速率用以初值为中心、分散度较大的对数正态（S=2）；`ucld.stdev` 用 `Normal(均值 0.4, sigma 0.3334)`（BEAUti 惯例），并按参数取值域截断；概率与频率参数用 `Uniform(下界, 上界)`；有符号参数用宽正态。离散状态节点（`IntegerParameter`、`BooleanParameter`）不需要先验。缺少这一步，`0` 到 `∞` 上的平坦先验就是非正常密度，分析也无法从 XML 重建。

**配置示例：**

```yaml
# Yule 模型
tree_prior:
  type: "yule"
  birth_rate: {value: 1.0, lower: 0.0}

# Birth-Death 模型
tree_prior:
  type: "birth_death"
  birth_rate: {value: 1.0, lower: 0.0}
  death_rate: {value: 0.5, lower: 0.0}
  sampling_rate: {value: 0.1, lower: 0.0}

# Coalescent 常量群体
tree_prior:
  type: "coalescent_constant"
  pop_size: {value: 10.0, lower: 0.0}

# Coalescent 指数增长
tree_prior:
  type: "coalescent_exponential"
  pop_size: {value: 10.0, lower: 0.0}
  growth_rate: {value: 0.1}

# Bayesian Skyline
tree_prior:
  type: "bayesian_skyline"

# CalibratedYule（需同时设置 calibration_method）
calibration_method: "calibrated_yule"
tree_prior:
  type: "calibrated_yule"
  birth_rate: {value: 1.0, lower: 0.0}
```

**在 `parameter_priors` 中为参数命名。** `parameter_priors:` 条目按生成器实际构建的状态节点匹配，别名下表：左列是配置中写的名字，右列是 XML 中出现的 id。

| 配置中写 | XML 中的 id |
|---|---|
| `birth_rate` | `birthRate`, `birthDiffRate` |
| `death_rate` | `relativeDeathRate` |
| `sample_probability`, `sampling_rate` | `sampleProbability`, `bdSkylineSamplingRate` |
| `pop_size` | `popSize` |
| `pop_sizes` | `bPopSizes`, `populationSizes` |
| `growth_rate` | `growthRate` |
| `clock_rate` | `clockRate` |
| `ucld_mean` / `ucld_stdev` | `ucld.mean` / `ucld.stdev` |
| `gamma_shape` / `proportion_invariant` | `gammaShape` / `proportionInvariant` |
| `kappa` | `kappa`（任意 partition 的 `…hky.kappa`） |

无法解析出目标状态节点的先验是**硬错误**，错误信息会列出当前被估计的状态节点，因为这样的先验会静默地永不参与求值。反之亦然：为 `estimate: false` 固定的参数指定先验同样是错误，固定参数时应删掉对应先验。

---

## 11. 校准分布

### 11.1 支持的分布类型

| 分布 | 配置名 | 参数 | 取值约束 | 典型用途 |
|---|---|---|---|---|
| Normal | `normal` | `mean`、`sigma` | `sigma > 0` | 约束良好的对称校准 |
| LogNormal | `lognormal` | `M`、`S`、`mean_in_real_space` | `S > 0`；仅当 `mean_in_real_space: true` 时要求 `M > 0` | 右偏校准（化石数据常用） |
| Uniform | `uniform` | `lower`、`upper` | 两者必填，且 `lower < upper` | 硬边界约束 |
| Exponential | `exponential` | `mean` | `mean > 0` | 最小年龄约束 |
| Gamma | `gamma` | `alpha`、`beta` | 两者均大于 0 | 灵活偏态分布 |
| Beta | `beta` | `alpha`、`beta` | 两者均大于 0 | 有界分布 |
| Laplace | `laplace` | `mu`、`scale` | `scale > 0` | 更重尾的对称分布 |
| InverseGamma | `inverse_gamma` | `alpha`、`beta` | 两者均大于 0 | 正实数上的右偏分布 |
| OneOnX | `one_on_x` | （无） | — | 尺度不变先验 |
| Poisson | `poisson` | `lambda` | `lambda > 0` | 计数数据 |
| ChiSquare | `chi_square` | `dof` | `dof > 0`（写入 BEAST2 的 `df` 输入） | 卡方分布 |

`uniform` 要过两道检查。参数本身必须满足 `lower < upper`（加上 `offset` 之后仍然成立），因为上下界颠倒的 `Uniform` 在整个支撑上密度为 0；`provenance.original_age_min` 也不得大于 `original_age_max`。

对 `lognormal` 而言，`M` 默认是**对数尺度**上的均值（中位数等于 `exp(M)`）。只有设置 `mean_in_real_space` 时 `M` 才是真实空间均值，此时 `M` 必须为正。

**配置示例：**

```yaml
# Normal 分布
distribution:
  type: "normal"
  parameters: {mean: 6.0, sigma: 0.5}

# LogNormal 分布（M 为 log 空间均值）
distribution:
  type: "lognormal"
  parameters: {M: 1.0, S: 0.5, mean_in_real_space: false}
  offset: 12.0    # 分布从 12.0 Ma 开始

# LogNormal 分布（mean_in_real_space=true 时 M 为真实空间均值）
distribution:
  type: "lognormal"
  parameters: {M: 15.0, S: 0.5, mean_in_real_space: true}

# Uniform 分布
distribution:
  type: "uniform"
  parameters: {lower: 20.0, upper: 40.0}

# Exponential 分布
distribution:
  type: "exponential"
  parameters: {mean: 5.0}

# Gamma 分布
distribution:
  type: "gamma"
  parameters: {alpha: 2.0, beta: 0.5}
```

### 11.2 Offset 支持

所有分布支持 `offset` 属性，用于平移分布原点。这在校准点有最小年龄约束时非常有用：

```yaml
distribution:
  type: "lognormal"
  parameters: {M: 1.0, S: 0.5}
  offset: 12.0    # 分布从 12.0 Ma 开始，即实际年龄 = offset + LogNormal(M, S)
```

### 11.3 超参数先验

校准分布参数本身可以设置先验分布（超参数先验）：

```yaml
calibrations:
  - name: "humanChimpMRCA"
    taxa: ["human", "chimp"]
    distribution:
      type: "normal"
      parameters: {mean: 6.0, sigma: 0.5}
      hyperpriors:
        mean:
          type: "uniform"
          parameters: {lower: 0.0, upper: 20.0}
        sigma:
          type: "uniform"
          parameters: {lower: 0.0, upper: 5.0}
```

每个超先验的键都必须是其所属分布的参数名：在 `normal`（参数为 `mean` 与 `sigma`）上写 `hyperpriors: {mu: …}`，解析器会拒绝。超先验自身的参数与任何校准分布一样接受值域检查。超先验会把该参数变成一个被估计的状态节点，这正是生成器同时为它配上算子与 `Prior` 的原因。

### 11.4 校准点来源元数据

每个校准点可附加来源元数据（`provenance`），用于记录化石来源和校准类型：

```yaml
provenance:
  source: "fossil"              # 来源类型：fossil / biogeographic / secondary / other
  reference: "DOI:10.1038/xxx"  # 参考文献 DOI 或 URL
  calibration_type: "soft"      # soft（软边界）/ hard（硬边界）
  original_age_min: 5.5         # 原始化石最小年龄
  original_age_max: 7.0         # 原始化石最大年龄
  notes: "化石描述文字"          # 自由文本备注
```

**`original_age_min/max` 现在会与所用先验的 95% 区间交叉核对**（区间由工具自带的 `compute_distribution_stats` 计算，`offset` 已计入）：

| 情形 | 结果 |
|---|---|
| `calibration_type: hard` 且先验 95% 区间完全落在记录年龄之外 | **错误** —— 硬边界校准必须真正包含它声称施加的年龄 |
| `calibration_type: soft` 且同样的不匹配 | **警告** —— 这是合理的：软边界本就刻意把质量放到化石下限之外（例如 2.5% 分位恰落在硬最小值处的对数正态），但应在方法学正文中说明 |
| `original_age_min > original_age_max` | **错误** |

```text
Calibration 'root': provenance records an original minimum age of 20.0 but the 95%
interval of the normal prior is [4.02, 5.98], i.e. entirely younger than the fossil.
A hard-bound calibration must contain the age it claims to enforce; either change
calibration_type to soft (with a soft-bound distribution) or move the prior.
```

### 11.5 校准点其他选项

```yaml
calibrations:
  - name: "rootCalibration"
    taxa: null                    # null 表示根节点校准
    monophyletic: true            # 强制单系约束（默认 true）
    use_originate: true           # 使用 stem 节点（originate）而非 MRCA
    tipsonly: false               # 是否仅末端校准（默认 false）
```

### 11.6 延伸阅读

- [CladeAge](https://github.com/BEAST2-Dev/cladeage)（BEAST2 插件）：基于 Birth-Death 速率与化石首现年龄，通过模拟生成校准密度。它是手工挑选 Normal/LogNormal 参数之外的另一种校准先验来源。
- [FossilCalibrations](https://fossilcalibrations.org)：经过同行评审的化石校准数据库；其 DOI 与年龄范围字段可直接对应上文的 `provenance` 元数据块。

---

## 12. 多 Partition 分析

### 12.1 多个序列文件

使用多个独立的序列文件作为不同 partition：

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
      linked_to: "gene1"    # 与 gene1 共享时钟模型
    tree: "shared"
```

### 12.2 密码子分区

使用 `FilteredAlignment` 实现密码子分区，基于同一序列文件的不同位点子集：

```yaml
alignments:
  - id: "alignment"
    file: "coding.fasta"
    data_type: "nucleotide"
  - id: "codon12"
    data_id: "alignment"       # 引用基础比对
    filter: "1::3,2::3"        # 第1+2位密码子
    data_type: "nucleotide"
  - id: "codon3"
    data_id: "alignment"
    filter: "3::3"             # 第3位密码子
    data_type: "nucleotide"

partitions:
  - id: "codon12"
    site_model:
      substitution_model: {type: "hky"}
      gamma_categories: 4
    clock_model:
      type: "ucln"
      ucld_mean: {value: 1.0}
      ucld_stdev: {value: 0.333}
    tree: "shared"
  - id: "codon3"
    site_model:
      substitution_model: {type: "gtr"}
      gamma_categories: 4
    clock_model:
      linked_to: "codon12"    # 共享 codon12 的时钟
    tree: "shared"
```

**Filter 表达式语法。** 表达式遵循 BEAST2 自身的 `FilteredAlignment` 语法，坐标为 **1-based**（从 1 开始计数）。表达式是一个逗号分隔列表：

| 形式 | 含义 |
|---|---|
| `N` | 单个位点 |
| `from-to` | 闭区间；空的 `from` 表示第 1 位，空的 `to` 表示最后一位 |
| `from-to\step` | 同一区间，每隔 `step` 个取一位 |
| `from:to:step` | 迭代器；`from:to:` 或 `from::step` 让空部分取默认值 |

于是 `1::3,2::3` 选出第 1、2、4、5、7、8 等位（密码子第 1、2 位），`3::3` 选出 3、6、9 等位。解析器会针对源比对逐项校验范围，以下情况都会报错，并指名违规项：无法解析的项、`step < 1`、超出比对长度的坐标、不选出任何位点的表达式。

**声明 filter 时 `data_id` 必填。** `data_id` 指向的源比对必须已先声明，且不得是过滤器比对自身：指向自身 id 会产出 `data="@self"` 并丢弃全部序列数据。上面的写法——源比对声明一次，每个被过滤分区用 `data_id` 指向它——正是 `examples/config_multi_partition.yaml` 与 `examples/config_advanced.yaml` 的用法，也正是让数据的未过滤副本恰好在 XML 中出现一次的机制。被过滤的比对若省略 `data_id`，解析器会报错（`… has a filter but no data_id: it would reference itself and lose all sequence data`）。

**两个分区不得覆盖同一批位点。** 解析器不会阻止同一份 FASTA 以两个 id 声明两次，并各建一个 `TreeLikelihood`。但由此得到的后验是两个相同似然之积，相当于把同一批 898 个位点计算两次，而 BEAST2 仍会照常运行。模型比较应拆成多次运行（每份配置一个候选模型，如 `examples/config_subst_jc69.yaml`、`config_subst_sym.yaml`、`config_subst_tim.yaml`、`config_subst_models.yaml`），分区分析则使用互不相交的 filter。

### 12.3 模型链接规则

- **`clock_model.linked_to: "partition_id"`**：该 partition 共享另一个 partition 的时钟模型。
- **`site_model.linked_to: "partition_id"`**：该 partition 共享另一个 partition 的位点模型（见 `config_advanced.yaml` 示例）。未指定时每个 partition 拥有独立的位点模型。
- `linked_to` 的值必须指向另一个已存在的 partition——自链接与未知 id 都是错误。

### 12.4 拓扑：只组装一棵连锁树

**`tree: "shared"` 是唯一合法取值。** 逐分区或独立拓扑尚未实现，因此解析阶段会直接拒绝该设置，而不是解析、存储之后再静默忽略：

```text
Partition 'codon3': tree: 'separate' requests a tree other than the shared one, but
Beast2Py assembles every partition on a single linked topology (per-partition /
unlinked trees are not implemented). Silently falling back to the shared tree would
turn an intended *BEAST-style multi-locus analysis into a fully linked one, so this
is an error. Remove the 'tree' key to analyse the shared topology, or run one
analysis per independent partition.
```

`tree: <某个 partition id>` 同样会遭到拒绝，错误信息一致。需要注意的后果是：本版本的多基因座分析假定单棵共享的物种树/基因树，每个 partition 的松弛时钟也索引同一棵树。独立基因树（`*BEAST`、StarDivergence 用例）在路线图上，本版本尚未实现。

---

## 13. 末端日期

末端日期（Tip Dates）用于末端采样时间各不相同的分析，如传染病传播研究。

```yaml
tip_dates:
  enabled: true
  trait_name: "date-forward"    # date-forward 或 date-backward，别无其他
  units: "year"                # year、month、day 或 coordinate
  dates:                        # 每个分类单元一项——覆盖率必须完整
    sample_2020: 2020
    sample_2019: 2019
    sample_2018: 2018
  # 替代/额外来源：两列 taxon,date 的 CSV
  # file: "tip_dates.csv"
```

**现已在解析阶段强制的规则**

- **完整覆盖。** `dates` 必须覆盖每个 partition 的每个分类单元。BEAST2 会把未定年的末端视为在时刻 0 采样，因此只给出部分映射时，缺失的支系就变成了同时刻采样的样本。解析器如今会准确报出缺失的分类单元（`tip_dates.dates covers 9 of 12 taxa; these are missing: …`），也会拒绝比对中不存在的分类单元。
- **`trait_name` ∈ {`date-forward`、`date-backward`}。** `date-forward` 表示从最早采样点向前计数（0 表示最古老样本），`date-backward` 从最晚采样点向后计数（0 表示最新样本，这也是 BEAUti 的 Tip dating 标签页产出的形式）。其他取值虽然也能产出合法 XML，但 BEAST2 不会把它们当作日期读取，分析就按同时刻采样运行。解析器会拒绝 `date_backward`（下划线写法），并在错误信息里列出合法名。
- **`units` ∈ {`year`, `month`, `day`, `coordinate`}。**
- **日期必须为数值。** CSV 中第二列若非数值，解析器会报错（表头行除外）。
- **保持单一时间轴。** `units` 和末端数值必须与校准共用同一时间轴。`units: year` 且日期为 0～20，旁边却是一个以 Ma 计的根校准，两者就相差六个数量级。要么全部用年表达（并相应缩放时钟速率），要么全部保持 Ma。

**会产出什么。** 启用末端日期后，生成器把 `TraitSet` 挂到树上，*并*为已定年的末端产出一个 `TipDatesRandomWalker` 算子，无需任何 `tipsonly` 校准。过去该算子受 `tipsonly` 校准门控：常见的 BEAUti 式“日期是性状而非校准”配置虽然产出 `TraitSet`，却没有算子移动它，运行因此返回一棵看起来已定年、实际未定年的树。

该算子只覆盖日期非零的末端。若全部末端都在同一时刻采样，就没有可移动的日期，也就不会产出 `TipDatesRandomWalker`——`examples/config_tip_dates.yaml` 正是这种情况，它的 12 个灵长类皆为现生样本。该行为已对 BEAST 2.7.8 验证：把其中一个末端改为早 1.5 Ma 的变体仍可解析并初始化。

---

## 14. 校准诊断

### 14.1 诊断概述

校准诊断是 Beast2Py 的核心创新功能，包含三个子模块：冲突检测器（`ConflictDetector`）、敏感性分析器（`SensitivityAnalyzer`）和可视化生成器（`VisualizationGenerator`）。

### 14.2 冲突检测

冲突检测器分析校准点对，检测三类问题：

**时间一致性冲突：** 对于嵌套 clade（分类集合 A ⊂ B），检查父 clade（B）的校准时间是否比子 clade（A）更老。比较 95% HPD 区间：父 clade 的 HPD 上界比子 clade 的 HPD 下界更年轻时，报告错误（`error`）；只有父 clade 的均值比子 clade 的均值更年轻、而 HPD 区间仍有重叠时，报告警告（`warning`）。

**分布重叠警告：** 计算校准分布之间的 1-Wasserstein 距离（Earth Mover's Distance）。对于非嵌套的校准点对，距离过小就报告警告，提示可能存在信息冗余。**默认判据是无量纲的**——Wasserstein 距离除以两个先验 95% 区间的平均宽度，低于 `relative_overlap_threshold`（默认 0.15）时警告。绝对判据（距离低于 `overlap_threshold`，默认 2.0 个时间单位）需显式开启，因为固定距离是尺度相关的：它随节点深度增长，因此在单一阈值下，浅节点几乎一定会触发，深节点则永不触发。如何切换与调参见[第 14.9 节](#149-diagnostics-块)。

**单系性冲突：** 检测有交集但不嵌套的分类集合。若两个校准都强制单系约束，这样的拓扑就不可能存在，检测器报告错误。

### 14.3 敏感性分析

敏感性分析器生成 N+1 个仅先验采样 XML（N 为校准点数量）：一个包含完整校准集的版本，以及 N 个留一法版本（逐个移除每个校准点）。用户用 BEAST2 运行这些短链 XML（100 万代即可），比较树高的后验分布，从而评估每个校准点的边际影响。（提示：仅先验采样本身是 BEAST2 的内置能力，对任意 XML 都可用 `beast -sampleFromPrior` 启用；敏感性分析器自动化的是 N+1 留一法设计。）这些文件落在报告 `diagnostics/` 目录下的 `sensitivity/` 子目录里；指定 `--sensitivity-xml-dir` 时，落在该目录的 `sensitivity/` 子目录里。

### 14.4 可视化

可视化生成器使用 matplotlib 创建所有校准分布的密度图，标注 95% HPD 区间（阴影区域）和均值（虚线）。支持 PNG、PDF 和 SVG 输出格式。

### 14.5 HTML 报告

诊断报告以自包含的 HTML 文件形式生成，包含分析概述、校准点摘要表、冲突检测结果表、分布可视化和敏感性分析说明。

### 14.6 CLI 运行诊断

```bash
beast2py diagnose \
    --config config.yaml \
    --report report.html \
    --sensitivity \
    --visualize \
    --sensitivity-xml-dir sensitivity_xmls/
```

### 14.7 API 运行诊断

```python
from beast2py.api import Beast2Py

b2p = Beast2Py()

# 完整诊断
report = b2p.diagnose(
    "config.yaml",
    output_dir="diagnostics/",
    run_sensitivity=True,
    run_visualization=True,
    report_path="report.html",
)

# 仅检测冲突
conflicts = b2p.detect_conflicts("config.yaml")

# 仅生成敏感性 XML
xmls = b2p.generate_sensitivity_xmls("config.yaml", output_dir="sensitivity/")

# 仅绘制分布图
b2p.plot_calibrations("config.yaml", output_file="calibrations.pdf")
```

### 14.8 冲突检测结果的解读

冲突检测返回 `Conflict` 对象列表，每个对象包含以下属性：`conflict_type`（冲突类型：temporal、distribution_overlap、monophyly）、`severity`（严重程度：error、warning、info）、`calibration_a` 和 `calibration_b`（涉及的校准点名称）、`description`（详细描述）和 `recommendation`（建议操作）。

- **error** 级别表示必须修复的问题（如时间不一致或单系性冲突）。
- **warning** 级别表示需要关注但非致命的问题（如分布重叠或均值不一致）。
- **info** 级别表示提示信息。

### 14.9 `diagnostics:` 块

`diagnostics:` 段可配置分布重叠检查的阈值。该段同样传给 `generate` 的闸门 2：

```yaml
diagnostics:
  overlap_measure: "relative"           # relative（默认）| absolute | both
  overlap_threshold: 2.0                # 绝对判据：以时间单位计的 Wasserstein 距离
  relative_overlap_threshold: 0.15      # 无量纲判据：W / 两先验 95% 区间平均宽度
```

- `overlap_measure` 决定启用哪些判据。`relative`（默认）即上文的尺度无关比值；`absolute` 恢复旧有的固定距离测试；`both` 在任一判据越界时告警。解析器会拒绝无法识别的值，并给出 `diagnostics.overlap_measure must be relative, absolute or both, got ...`。
- `overlap_threshold`（默认 `2.0`）与 `relative_overlap_threshold`（默认 `0.15`）是两个阈值，只有 `overlap_measure` 点名的判据才会生效。

默认改用无量纲判据后，固定的 `2.0` 在同一棵树内节点年龄悬殊时就不再含义一致。如今警告只在两个先验*相对于自身宽度*几乎可以互换时触发，这才是冗余校准信息的尺度无关含义。

---

## 15. 可重复性工具

### 15.1 分析指纹

分析指纹是基于*科学性*配置 SHA-256 哈希生成的确定性标识符，采用不含日期的格式 `B2P-{hash12}-{version}`，例如 `B2P-ba18c26ed61f-0.1.0`。它**不含日期与时区**，因此同一配置在同一比对数据上每次重跑都逐字节复现。

以下变化会改变指纹（`models.py::to_dict` 把这些都纳入哈希）：序列内容、树设置（`Partition.tree`）、冠群/茎群（`use_originate`）或 `tipsonly` 标志、超先验、`parameter_priors` 条目、末端日期、算子权重。只改输出文件名**不会**改变指纹。

除标识符外，XML 注释还记录一条独立的各比对内容摘要：`<!-- Analysis Fingerprint: B2P-{hash12}-{version} | Data: {hash16} -->`，其中 `Data:` 是对各比对序列摘要再取的前 16 位。指纹嵌入为 XML 注释，也可写入 JSON 侧车文件。墙钟时间 `generation_time` **只**出现在侧车文件中，绝不出现在标识符或 XML 里。实测 JSON 内容见[第 5.7 节](#57-fingerprint--生成分析指纹)。

```bash
# CLI 方式
beast2py fingerprint --config config.yaml --output fingerprint.json --verbose

# 或在生成时同时生成
beast2py generate --config config.yaml --output output.xml --fingerprint
```

```python
# API 方式
fp = b2p.generate_fingerprint("config.yaml", output="fingerprint.json")
print(f"分析指纹: {fp}")
```

### 15.2 方法学描述

自动生成 LaTeX 格式的方法学段落，描述替换模型（含引用）、时钟模型、树先验、校准点（含分布参数和来源）、MCMC 设置和分析指纹。可直接插入论文手稿。

```bash
beast2py methods --config config.yaml --output methods.tex
```

```python
methods_text = b2p.generate_methods("config.yaml", output="methods.tex")
```

### 15.3 管道文件生成

生成 Snakemake 或 Nextflow 管道文件，包含 XML 生成、BEAST2 执行、诊断和方法学描述的完整规则或进程。

```bash
# Snakemake
beast2py pipeline --config config.yaml --type snakemake --output-dir pipeline/

# Nextflow
beast2py pipeline --config config.yaml --type nextflow --output-dir pipeline/

# 两者都生成
beast2py pipeline --config config.yaml --type both --output-dir pipeline/
```

```python
# API 方式
path = b2p.generate_pipeline("config.yaml", pipeline_type="snakemake", output_dir="pipeline/")
sm_path, nf_path = b2p.generate_pipeline("config.yaml", pipeline_type="both", output_dir="pipeline/")
```

生成的 Snakemake 文件为 `Snakefile`，Nextflow 文件为 `main.nf`。

---

## 16. 验证

### 16.1 结构验证

结构验证（`validator.py::XMLValidator`）在一次调用里跑两层。

**结构层**检查以下各项：

- XML 格式正确，根元素为 `<beast>`
- 必需元素（`<data>` 与 `<run>`）与 `<state>` 元素存在
- id/idref/`@ref` 引用一致，无重复 ID
- 后验引用与日志器存在
- 命名空间属性含 `beast.base`

**语义层**（`validator.py::_check_semantic`）负责让分析*正确*，而不只是可解析：

- 每个被估计状态节点都带有 `Prior`。无界支撑上的平坦先验是非正常密度，也无法从 XML 重建；有界支撑的离散 `IntegerParameter`/`BooleanParameter` 节点豁免。
- 每个被估计状态节点都配有算子，且没有算子挂在固定参数上。
- 状态节点上下界不反向（`lower < upper`）。
- 初值落在自身界内，且与声明维度匹配。
- 每个 `<data>` 块的序列非空、等长，分类单元名唯一。

`storeEvery` 关闭续算时，给出警告而非错误。

```bash
beast2py validate --xml output.xml
```

```python
result = b2p.validate_xml("output.xml")
print(f"有效: {result.is_valid}")
print(f"错误: {result.errors}")
print(f"警告: {result.warnings}")
```

### 16.2 BEAST2 原生验证

`--beast2` 使用 BEAST2 的 `XMLParser` 做更严格的验证。Beast2Py 随包发布一个无头验证工具（`beast2py/tools/Beast2Validator.java`），无需 JavaFX 环境即可运行。它解析文件**并对每个对象调用 `initAndValidate()`**，因此检查的是模型能否**完成初始化**，不只是能否解析。拼错的可选属性、永不触发的算子、非法的约束都会在此失败，而不是照样打印 `VALID`。

该 Java 助手及其 `launcher.jar`、`classes/` 位于包内，作为 `package-data` 随 wheel 发布。`validator.py` 会在安装包旁、项目根、`PATH` 上或 `--beast2-path` 指定处定位 `beast2_validate.sh`。

```bash
beast2py validate --xml output.xml --beast2 --beast2-path /path/to/beast
```

```python
result = b2p.validate_xml("output.xml", beast2_validate=True, beast2_path="/path/to/beast")
```

这条路径与 `generate` 的闸门 3 完全相同，需要安装 BEAST2 与一个 **JDK 17 或更新版本**，`--beast2-path` 可指向特定安装或脚本。随包脚本在[第 2.6 节](#26-安装-beast2可选用于-beast2-验证闸门)所列的候选中搜索该 JDK，而不是只信 `PATH`：`PATH` 上一个更旧的 `java` 不该遮蔽可用的 `openjdk@17`。所有候选都不合格时，脚本给出清晰的信息并以退出码 2 放弃。

`beast2_validation_status` 区分三种结果：`passed`、`failed` 与 `unavailable`。缺少 BEAST2 安装报告为*不可用*，绝不会报告为通过。由于这项检查是显式要求的，生成会以退出码 2 停止，除非加 `--allow-unvalidated`。

### 16.3 验证状态

全部 **19** 个示例配置文件均通过三道闸门：结构检查、校准冲突检测，以及 BEAST2 `parseFile` 与 `initAndValidate` 检查。第三道闸门使用无头模式下的原生 BEAST2 **v2.7.8** 解析器，无需 JavaFX。其中两份依赖附加包的文件——`config_bd_skyline.yaml` 与 `config_nested_sampling.yaml`——已对照 `~/.beast/2.7/` 下安装的包核实，即 **BDSKY 1.5.1** 与 **NS 1.2.0**。

测试套件共收集 **359** 个测试：**338** 个单元/语义/发布完整性测试与 **21** 个 BEAST2 集成测试。集成测试对全部 19 个示例配置、`quick` 路径以及仓库内已提交的 `output_basic.xml` 生成 XML，并用真实的 BEAST 2.7.8 解析器与模型初始化逐一检查。当环境缺少 BEAST2 或 JDK 17 时，这些测试直接跳过，不算失败。全部通过。

---

## 17. MCMC 设置

### 17.1 标准 MCMC

```yaml
mcmc:
  type: "standard"
  chain_length: 10000000      # 链长度（总迭代次数）
  pre_burnin: 0               # 预烧入迭代数（必须 < chain_length）
  store_every: 1000           # 检查点间隔，按*已记录样本数*计（见下）
  sample_from_prior: false    # 仅从先验采样（用于敏感性分析）
  seed: 1                     # 供 `beast -seed` 使用的随机种子（**不写入 XML**，见下）
```

**`store_every` 计量的是已记录的样本数，不是代数。** BEAST2 的 `MCMC` 判断 `(sampleNr + 1) % storeEvery == 0`（`MCMC.java:490`）。因此在默认的 `logEvery` 下，`store_every` 是两个 `.xml.state` 检查点之间*保存的样本*数，也就是中断后可续算的位置个数。默认值**随链长伸缩**，即 `max(1, chain_length // 10000)`：链长 1,000 万时为 `1000`，100 万时为 `100`。这样任何一次运行都能得到大约十个检查点，而不是只有最末的一个。显式写 `-1` 仍然合法，但会关闭检查点，此时验证器给出 `storeEvery='-1' ... cannot be resumed` 警告。

**生成器刻意不把 `mcmc.seed` 写进 XML。** BEAST2 的 `MCMC` 没有 `seed` *输入*，写出一个会让模型无法解析（`has no input with name seed`）。随机数种子是命令行参数（`beast -seed N`）。Beast2Py 把 `mcmc.seed` 保留在配置里，供 CLI 与启动器消费。日志文件名中仍可携带 `$(seed)` 变量，由 BEAST2 在运行时展开（见[第 18 节](#18-日志器配置)）。

### 17.2 Nested Sampling MCMC

Nested Sampling 用于模型比较和边缘似然估计，需要 BEAST2 NS add-on 包。

```yaml
mcmc:
  type: "nested_sampling"
  chain_length: 10000000      # 每个子链的长度
  particle_count: 1           # 粒子数
  sub_chain_length: 10000     # 子链长度
```

### 17.3 仅先验采样

将 `sample_from_prior` 设为 `true` 可生成仅从先验采样的 XML（不使用序列数据）。这通常由敏感性分析器自动设置，也可手动使用：

```yaml
mcmc:
  type: "standard"
  chain_length: 1000000       # 短链即可
  sample_from_prior: true
```

**提示：** 对任意既有 BEAST2 XML，也可以不修改文件而直接通过 BEAST2 命令行实现同样效果：`beast -sampleFromPrior output.xml`。Beast2Py 的增值在于自动生成完整的 N+1 个留一法敏感性 XML（见第 14.3 节）。

---

## 18. 日志器配置

日志器控制 BEAST2 运行时输出文件的频率和格式。

```yaml
loggers:
  trace_log:                    # 参数轨迹日志（.log 文件）
    file_name: "output.$(seed).log"
    log_every: 1000             # 每 1000 代记录一次
  tree_log:                     # 树日志（.trees 文件）
    file_name: "output.$(seed).trees"
    log_every: 1000
  screen_log:                   # 屏幕输出
    log_every: 10000            # 每 10000 代输出一次
```

**`$(seed)` 变量：** BEAST2 会把文件名中的 `$(seed)` 替换为随机种子值，确保多次运行不会互相覆盖。

---

## 19. 初始化策略

Beast2Py 支持三种初始树生成策略：

```yaml
initialization:
  tree_type: "random"     # 随机树（默认）
  # tree_type: "upgma"    # UPGMA 聚类树
  # tree_type: "newick"   # 用户提供的 Newick 树
  newick_file: null        # newick 类型时指定树文件路径
```

- **`random`**：生成随机拓扑树，适用于大多数分析。
- **`upgma`**：基于序列距离生成 UPGMA 聚类树，提供更好的起始点。
- **`newick`**：使用用户提供的 Newick 格式树文件，适用于需要特定起始拓扑的分析。

---

## 20. 操作器权重

用户可以选择性地覆盖默认的 MCMC 操作器权重。树操作器的默认权重为：`tree_scaler` 3.0、`uniform` 30.0、`subtree_slide` 15.0、`exchange_narrow` 15.0、`exchange_wide` 3.0、`wilson_balding` 1.0、`up_down` 3.0。

```yaml
operator_weights:
  tree_scaler: 3.0        # 树缩放操作（默认 3.0）
  uniform: 30.0           # 均匀操作（默认 30.0）
  subtree_slide: 15.0     # 子树滑动（默认 15.0）
  exchange_narrow: 15.0   # 窄交换（默认 15.0）
  exchange_wide: 3.0      # 宽交换（默认 3.0）
  wilson_balding: 1.0     # Wilson-Balding 操作（默认 1.0）
  # ... 其他操作器
```

该功能面向高级用户，通常无需修改默认权重。

---

## 21. 示例配置文件

`examples/` 目录包含 **19** 个完整示例配置文件。如实统计，它们用到了**全部 4 种时钟模型**（strict、UCLN、UCE、RLC）、**8 种树先验中的 7 种**（除 `coalescent_constant` 外全部），以及**13 种替换模型中的 9 种**——7 种核苷酸模型（JC69、HKY、TN93、GTR、SYM、TIM、TVM）加两种氨基酸矩阵（WAG、JTT）。未被示例覆盖的四种是氨基酸的 Dayhoff、BLOSUM62、CPREV 与 MTREV。它们**并不**把“所有模型塞进一个文件”：替换模型的比较拆成 `config_subst_jc69/sym/tim/models.yaml` 四份配置，一次运行一个模型。原因是在同一棵共享拓扑上放多个似然，得到的是后验的自乘，而不是模型比较（见第 12.2 节）。

| 文件 | 场景 | 树先验 | 替换模型 | 时钟模型 | 校准数 | 特殊功能 |
|---|---|---|---|---|---|---|
| `config_basic.yaml` | 基础分析 | Yule | HKY+Γ4 | Strict | 1 | 最简配置；显式 `parameter_priors` + 来源元数据；`output_basic.*` 的来源 |
| `config_advanced.yaml` | 高级多校准 | Birth-Death | GTR+Γ4+I / 链接位点模型 | UCLN + 链接 | 3 | 密码子分区 + 位点*与*时钟链接 + 超先验 + `diagnostics` 段 + stem 校准 |
| `config_subst_models.yaml` | 替换模型（TVM） | Yule | TVM+Γ4 | Strict | 1 | 一次运行一个模型：4 个颠换率 + 1 个共享转换率 |
| `config_subst_jc69.yaml` | 替换模型（JC69） | Yule | JC69+Γ4 | Strict | 1 | 等速率/等频率基线模型 |
| `config_subst_sym.yaml` | 替换模型（SYM） | Birth-Death | SYM+Γ4 | UCLN | 2 | 六个自由交换率 + `frequencies.mode: uniform` |
| `config_subst_tim.yaml` | 替换模型（TIM） | Yule | TIM+Γ4 | Strict | 1 | rateAG/rateCT + 两个颠换率，`frequencies.mode: empirical` |
| `config_aa_model.yaml` | 氨基酸模型 | Birth-Death | WAG+Γ4 | UCLN | 2 | 氨基酸数据 + Laplace/LogNormal 校准 |
| `config_aa_jtt.yaml` | 氨基酸模型（JTT） | Yule | JTT+Γ4+I | Strict | 2 | 蛋白数据上的 Gamma+I |
| `config_tip_dates.yaml` | 末端日期 | Birth-Death | HKY | UCLN | 1 | `TraitSet` 末端日期（12 个分类单元全部给日期） |
| `config_calibrated_yule.yaml` | CalibratedYule | CalibratedYule | HKY+Γ4 | Strict | 3 | `calibration_method` 切换 + stem 校准 |
| `config_multi_partition.yaml` | 多 partition | Birth-Death | HKY+Γ4 / GTR+Γ4 | UCLN + 链接 | 2 | `FilteredAlignment` 密码子分区 + 时钟链接，单棵连锁拓扑 |
| `config_relaxed_clock.yaml` | 松弛对数正态时钟 | Birth-Death | GTR+Γ4 | UCLN | 2 | 有界的 `ucld_mean` / `ucld_stdev` |
| `config_relaxed_exponential.yaml` | 松弛指数时钟 | Yule | GTR+Γ4 | UCE | 2 | UCE 松弛时钟 |
| `config_rlc.yaml` | 随机局部时钟 | Yule | HKY+Γ4 | RLC | 1 | 指示变量与速率状态节点 |
| `config_coalescent_exponential.yaml` | 指数增长溯祖 | Coalescent（指数） | HKY+Γ4 | Strict | 1 | 有符号的 `growth_rate` |
| `config_skyline.yaml` | 贝叶斯天际线 | Bayesian Skyline | TN93+Γ4 | Strict | 1 | 天际线群体历史 |
| `config_ebsp.yaml` | EBSP | EBSP | HKY+Γ4 | Strict | 1 | 扩展贝叶斯天际线图（BEAST 2.7 核心包内置） |
| `config_bd_skyline.yaml` | BD Skyline Serial | BD Skyline Serial | HKY+Γ4 | Strict | 1 | 需 bdsky 插件包（已对照 BDSKY 1.5.1 核实） |
| `config_nested_sampling.yaml` | Nested Sampling | Yule | HKY+Γ4 | Strict | 1 | NS MCMC（需 NS 插件包，已对照 NS 1.2.0 核实） |

表中的“校准数”是每个文件里 `calibrations:` 条目的实际个数，模型列是各分区实际声明的内容（`/` 分隔不同分区）。

此外，`calibrations.yaml` 是供 `quick` 子命令使用的独立校准点文件。示例数据包括 `primates.fasta`（12 个分类单元 × 898 bp）和 `aminoacid.fasta`（10 个分类单元 × 234 aa），示例输出制品包括 `output_basic.xml`、`output_basic.fingerprint.json` 和 `output_basic.methods.tex`。

**运行示例：**

```bash
# 运行基础示例
beast2py generate -c examples/config_basic.yaml -o output_basic.xml -v

# 运行高级示例（含诊断）
beast2py generate -c examples/config_advanced.yaml -o output_advanced.xml --diagnose --fingerprint --methods -v

# 运行多 partition 示例
beast2py generate -c examples/config_multi_partition.yaml -o output_mpart.xml -v

# 快速模式示例
beast2py quick -a examples/primates.fasta -o output_quick.xml \
    --tree-prior yule --subst-model hky --clock-model strict \
    --calibration-yaml examples/calibrations.yaml
```

```python
# API 运行示例
from beast2py.api import Beast2Py

b2p = Beast2Py(base_dir="examples/")

# 基础示例
xml = b2p.generate_xml("config_basic.yaml", output="output_basic.xml")

# 高级示例（含所有功能）
result = b2p.generate(
    "config_advanced.yaml",
    output="output_advanced.xml",
    diagnose=True,
    fingerprint=True,
    methods=True,
)
```

---

## 22. 故障排除

### 22.1 严格校验：被拒绝的输入及其错误信息

Beast2Py 是**失败即关闭**的：生成器无法兑现某个配置时，`ConfigError` 会在写出任何文件*之前*抛出（CLI：`[ERROR] Configuration error: …`，退出码 1；Python API：抛出异常）。下面这些是当前代码实际打印的字符串。

**Schema 与未知的键。** 解析器把每个段都对照一份键白名单校验（`config.py::KNOWN_KEYS`）。不属于该段 schema 的键会报错，而不是静默地填入默认值：

```text
mcmc: unknown key(s) chainlength_typo. Allowed keys for this section: chain_length,
particle_count, pre_burnin, sample_from_prior, seed, store_every, sub_chain_length,
type. Beast2Py refuses unknown keys rather than silently applying a default, because
a typo in a model or MCMC setting changes the analysis without changing the XML's
validity.
```

**布尔值。** 解析器只接受真正的 YAML 布尔值，或 `true/false/yes/no/on/off/1/0` 之一，并拒绝任何含混的值（带引号的 `"false"` 仍然接受，含义为*假*）：

```text
'alignment.hky.kappa': 'maybe' is not a boolean; use true/false (unquoted in YAML) or
one of ['0', '1', 'false', 'no', 'off', 'on', 'true', 'yes']
```

**没有类别数的 Gamma 形状**（与 `SiteModel.java:64` 一致，即 `Ignored if gammaCategoryCount 1 or less`）：

```text
Partition 'alignment': gamma_shape was supplied but gamma_categories is 1. SiteModel.java
documents the shape input as 'Ignored if gammaCategoryCount 1 or less', so an estimated
shape would be a pure prior random walk that contributes nothing to the likelihood while
still appearing as a column in the trace.
```

**固定时钟速率且没有校准**（绝对时间不可识别）：

```text
The clock rate is fixed and no calibration is supplied, so absolute time is not
identifiable: node heights can only be estimated in expected-substitutions-per-site
units. Either estimate the clock rate or add at least one absolute (e.g. fossil)
calibration.
```

**独立树与逐分区的树**（`tree: separate` 或一个 partition id）：

```text
Partition 'alignment': tree: 'separate' requests a tree other than the shared one, but
Beast2Py assembles every partition on a single linked topology (per-partition /
unlinked trees are not implemented). ... Remove the 'tree' key to analyse the shared
topology, or run one analysis per independent partition.
```

**末端日期覆盖度**（每个分类单元都必须有日期，未知或缺失的分类单元会报错）：

```text
tip_dates.dates covers 2 of 12 taxa; these are missing: Gorilla, Homo_sapiens, ...
BEAST2 treats an undated tip as sampled at time zero, which silently turns those
lineages into contemporaneous samples. Supply a date for every tip, or disable tip_dates.
```

**MCMC 长度规则**（`chain_length > 0`、`pre_burnin >= 0` 且 `< chain_length`、`store_every` 为正数或 `-1`）：

```text
mcmc.pre_burnin (5000) must be smaller than chain_length (100); otherwise the chain
never samples.
```

**Filter 语法与范围**（被过滤的比对必须有 `data_id`，1-based 的区间会对照源比对长度校验）：

```text
Alignment 'codon3': 'filter' requires 'data_id' naming the source alignment; without it
the filter would reference itself
```

**频率向量必须合法**（长度与 `dimension` 匹配、逐项大于 0、在 2% 容差内和为 1，随后归一化）：

```text
alignment.frequencies: frequencies sum to 1.2, not 1.0. Supply a proper probability
vector (a rounding tolerance of 0.02 is accepted and the vector is then renormalised).
```

**反向的校准边界：**

```text
Calibration 'root': uniform lower (40.0) must be strictly below upper (20.0); as given,
the density is zero on the whole support.
```

**生灭过程的定义域**（`death_rate < birth_rate`；`death >= birth` 请改用 BDSKY）：

```text
tree_prior: death_rate (0.8) must be smaller than birth_rate (0.5) for the
BirthDeathGernhard08Model, whose parameters are birth - death (>= 0) and death / birth
(< 1). A super-critical process (death >= birth) needs the BDSKY add-on
(tree_prior.type: bd_skyline_serial), which parameterises birth, death and sampling
rates separately.
```

**无法解析的 `parameter_priors` 名称**（匹配不到任何被估计节点的先验会报错，而不是静默跳过）：

```text
parameter_priors: 'totally_bogus_name' does not match any estimated parameter of this
analysis. Estimated state nodes are: alignment.clockRate, alignment.freqParameter,
alignment.hky.kappa, birthRate. A prior whose target does not exist is silently never
evaluated.
```

### 22.2 常见错误

**`Configuration file not found`**

确保 YAML 文件路径正确。使用绝对路径或相对于配置文件所在目录的相对路径。序列文件路径相对于 YAML 配置文件所在目录解析。

**`Invalid YAML in <file>: …`**

文件根本无法解析：流式序列未闭合、用制表符做缩进、键名重复。消息会点出文件名，并保留解析器
自己的行列号，因此要改的位置直接写在消息里，而不是留给栈回溯去暴露。所有读取配置的命令都按
这一形式报错，`quick --calibration-yaml` 也一样。

**`Calibration file not found: …`**

`quick --calibration-yaml` 需要的是一份校准点 **列表** 形式的 YAML，而不是完整的分析配置。
若文件顶层是映射而非列表，会以同样形式被拒绝（`Calibration YAML must be a list of calibration points`）。

**`Calibration taxon not found in any alignment`**

检查校准中的分类单元名称是否与序列头完全匹配。FASTA 序列头取 `>` 之后的第一个词作为分类单元名，且大小写敏感。

**`Partition ID not matching alignment`**

确保 partition ID 与比对 ID 匹配。如果使用 `FilteredAlignment`（密码子分区），partition ID 应与过滤后比对的 ID 匹配，而不是基础比对的 ID。

**`Duplicate partition ID`**

每个 partition 必须有唯一 ID。检查 YAML 配置中是否有重复的 partition ID。

**`Invalid substitution model`、`Invalid clock model`、`Invalid tree prior type`**

使用 `beast2py list-models` 查看支持的模型类型，确保配置中使用的名称正确。

**`Namespace does not contain 'beast.base'`**

此警告表示生成的 XML 使用 BEAST2 v2.6.x 命名空间格式。Beast2Py 生成 v2.7.x 格式。如果使用 BEAST2 v2.6.x，需要升级到 v2.7.0 或更新版本。

**`BEAST2 executable not found`**

确保 BEAST2 已安装且 `beast` 命令在 `PATH` 上，或通过 `--beast2-path` 指定完整路径。

**`Could not find package: bdsky / NS`（插件包缺失）**

BD Skyline Serial 与 Nested Sampling 示例需要 [bdsky](https://github.com/BEAST2-Dev/bdsky) 和 [NS](https://github.com/BEAST2-Dev/nested-sampling) 插件包。可通过 BEAST2 Package Manager 图形界面安装，或使用命令行（`launcher.jar` 为 BEAST2 自带 `lib` 目录中的同名文件）：

```bash
java -cp launcher.jar beast.util.PackageManager -add bdsky
java -cp launcher.jar beast.util.PackageManager -add NS
```

安装后用 Package Manager 或 `beast -version` 确认插件版本；R 用户也可使用 [mauricer](https://github.com/ropensci/mauricer) 包的等价辅助函数（`is_beast2_pkg_installed()`、`get_beast2_pkg_names()`）。

### 22.3 诊断中的常见问题

**时间一致性错误**

父 clade 的校准时间比子 clade 更年轻。检查嵌套校准点的分布参数，确保父 clade 的 HPD 区间覆盖的时间范围比子 clade 更老。

**分布重叠警告**

按默认的（相对）判据，检测器计算两个非嵌套校准之间的 1-Wasserstein 距离，除以它们 95% 区间的平均宽度，该比值低于 `relative_overlap_threshold`（默认 0.15）时报告警告。也就是说，两个先验*相对于自身宽度*几乎可以互换。检查这些校准是否提供了冗余信息，能否合并。若要切换回（或叠加）绝对的 `overlap_threshold` 判据，见[第 14.9 节](#149-diagnostics-块)。

**单系性冲突**

两个有交集但不嵌套的分类集合同时要求单系。确保校准的分类单元集合要么完全不相交，要么完全嵌套；或放宽其中一个校准的单系约束（设置 `monophyletic: false`）。

### 22.4 API 使用问题

**`ImportError: No module named 'beast2py'`**

确保已通过 `pip install -e .` 安装了包，且当前工作目录或 Python 路径包含项目目录。

**API 调用时的 `FileNotFoundError`**

如果配置中使用了相对路径的序列文件，确保 `base_dir` 参数指向正确的目录：`Beast2Py(base_dir="examples/")`。

### 22.5 获取帮助

```bash
# 一般帮助
beast2py --help

# 命令特定帮助
beast2py generate --help
beast2py diagnose --help
beast2py quick --help
beast2py validate --help
beast2py methods --help
beast2py pipeline --help
beast2py fingerprint --help
beast2py list-models --help
```

```python
# Python 帮助
import beast2py
help(beast2py.api)
help(beast2py.api.Beast2Py)
```

### 22.6 报告问题

请在 https://github.com/ZengZichao/Beast2Py/issues 报告错误和功能请求。提交问题时请附上：Beast2Py 版本（`beast2py --version`）、Python 版本、操作系统、配置文件内容（如适用）、完整错误信息。

---

## 23. 附录：模型完整列表

### 23.1 替换模型（13 种）

**核苷酸模型（7 种）：** `jc69`、`hky`、`tn93`、`gtr`、`sym`、`tim`、`tvm`

**氨基酸模型（6 种）：** `wag`、`jtt`、`dayhoff`、`blosum62`、`cprev`、`mtrev`

### 23.2 时钟模型（4 种）

`strict`、`ucln`、`uce`、`rlc`

### 23.3 树先验（8 种）

`yule`、`calibrated_yule`、`birth_death`、`coalescent_constant`、`coalescent_exponential`、`bayesian_skyline`、`ebsp`、`bd_skyline_serial`

### 23.4 校准分布（11 种）

`normal`、`lognormal`、`uniform`、`exponential`、`gamma`、`beta`、`laplace`、`inverse_gamma`、`one_on_x`、`poisson`、`chi_square`

### 23.5 数据类型（2 种）

`nucleotide`、`aminoacid`

### 23.6 MCMC 类型（2 种）

`standard`、`nested_sampling`

### 23.7 初始化策略（3 种）

`random`、`upgma`、`newick`

### 23.8 校准方法（2 种）

`mrca_prior`、`calibrated_yule`
