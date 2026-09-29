# Beast2Py

**一个用于可重复分歧时间估计的 Python 框架，具备自动化校准先验指定、验证和诊断功能。**

[![许可证: MIT AND LGPL-2.1-only](https://img.shields.io/badge/许可证-MIT%20AND%20LGPL--2.1--only-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![BEAST2 2.7.8](https://img.shields.io/badge/BEAST2-2.7.8-orange.svg)](https://www.beast2.org/)

---

## 概述

Beast2Py 是一个 Python 框架，用于**规范化指定、验证和诊断 BEAST2 分歧时间估计中的节点校准先验**，并生成可重复的 BEAST2 分析配置。该工具同时提供命令行接口（CLI）和 Python 应用程序编程接口（API），使研究人员能够从人类可读的 YAML 配置文件生成完整的 BEAST2 XML 文件，自动检测校准冲突，执行先验敏感性分析，并生成可重复性制品。

### 核心创新

> 功能对比基线：BEAUti2 v2.7.8（BEAST.app 包）与 beautier v2.6.15。

| 功能 | 描述 |
|------|------|
| **校准冲突检测** | 自动检测多个校准点之间的时间不一致性、分布重叠和单系性冲突 |
| **先验敏感性分析** | 自动生成“仅先验采样”XML，评估校准先验对后验的边际影响 |
| **校准来源元数据** | 每个校准点记录化石来源、参考文献 DOI、校准类型（硬/软边界） |
| **可重复性框架** | YAML 配置版本化 + 分析指纹 + LaTeX 方法学描述 + Snakemake/Nextflow 管道生成 |

### 三层功能集

- **第一层（beautier 对齐）：** R 包 beautier 支持的全部功能
- **第二层（beautier 不支持）：** beautier “Missing features/unsupported” 清单中计入二级子项后的全部 13 项（该清单有 10 个顶级条目；“Clock models” 按其 2 个子项计、“Tree priors” 按其 3 个子项计，故 8 + 2 + 3 = 13）：分布 offset、两条以上 DNA 比对、两套以上位点/时钟/树模型、两个以上 MRCA 先验、共享位点/时钟/树模型、氨基酸比对、超参数、松弛指数时钟、随机局部时钟、Calibrated Yule、Coalescent EBSP、Birth Death Skyline Serial、初始化策略。在这 13 项之外，Beast2Py 另加 SYM、TIM、TVM 核苷酸替换模型
- **第三层（独有整合）：** 7 项诊断与可重复性能力在配置阶段的独有整合。个别能力在其他工具中以孤立形式存在（如 BEAUti2 2.7+ 内置交互式方法学段落查看器与先验分布预览、`beast -sampleFromPrior` 可对任意 XML 做先验采样），但作为一体化、可脚本化的配置阶段工作流尚无先例

## 支持的功能

- **序列数据：** FASTA / NEXUS 导入（核苷酸和氨基酸）
- **替换模型：** JC69、HKY、TN93、GTR、SYM、TIM、TVM，另有氨基酸模型（WAG、JTT、Dayhoff、Blosum62、CPREV、MTREV）
- **位点模型：** Gamma 速率异质性、不变位点、Gamma+I
- **时钟模型：** 严格时钟、UC 松弛对数正态（UCLN）、UC 松弛指数（UCE）、随机局部时钟（RLC）
- **树先验：** Yule、CalibratedYule、Birth-Death、Coalescent（常量/指数）、Bayesian Skyline、EBSP；BD Skyline Serial（需要 [bdsky](https://github.com/BEAST2-Dev/bdsky) 插件包）
- **校准分布：** Normal、LogNormal、Uniform、Exponential、Gamma、Beta、Laplace、InverseGamma、OneOnX、Poisson、ChiSquare
- **多 Partition：** 多条比对（或同一比对的 `FilteredAlignment` 密码子位点子集），位点模型与时钟可独立或共享（`linked_to`）。所有 partition 都装配在**同一棵连锁拓扑**上：解析阶段会直接拒绝 `tree: separate` 与任何 partition 专属树 id，并给出解释，因为不连锁基因树（`*BEAST`、StarDivergence 一类的分析）尚未实现
- **多 MRCA 先验：** 无限校准点
- **末端日期：** 用 `TraitSet` 支持时序采样数据（`TipDatesRandomWalker` 在启用末端日期时即产出，不再依赖 `tipsonly` 校准）
- **MCMC：** 标准 MCMC + Nested Sampling MCMC（嵌套抽样需要 [NS](https://github.com/BEAST2-Dev/nested-sampling) 插件包）
- **初始化：** RandomTree、UPGMA、User Newick
- **超参数先验：** 为校准分布的参数设置先验
- **严格配置校验：** 解析阶段会拒绝未知键、加引号的布尔值、不等长或重名比对、标错的数据类型、定义域外的先验以及不可辨识的时间轴，而不是静默地采用默认值

## 安装

### 前置要求

- Python 3.10 或更高版本
- BEAST2 v2.7.8（已核实的版本；生成的 XML 面向 BEAST 2.7 的 `beast.base.*` 命名空间，其他 2.7.x 预期可用但未实测）。BEAST 2.7.x 的 class 文件要求 **JDK 17 或更新版本**
- JDK 17 或更新版本，仅在使用可选的 BEAST2 验证闸门（`--beast2-validate` / `validate --beast2`）时需要；其余功能为纯 Python

### 从源码安装

```bash
git clone https://github.com/ZengZichao/Beast2Py.git
cd Beast2Py
pip install -e .
```

### 依赖

```
pyyaml >= 6.0
biopython >= 1.80
numpy >= 1.20
scipy >= 1.8
matplotlib >= 3.5
defusedxml >= 0.7
```

## 快速开始

### 命令行使用

以下每条命令都在仓库根目录实测通过。`pip install -e .` 之后可直接用 `beast2py …`；
未安装时把 `beast2py` 换成 `python3 -m beast2py.main` 即可（行为一致）。

```bash
# 从 YAML 配置文件生成 BEAST2 XML
beast2py generate --config examples/config_basic.yaml --output output_basic.xml -v

# 运行校准诊断（冲突检测 + 敏感性分析 + 可视化）
beast2py diagnose --config examples/config_basic.yaml --report report.html --sensitivity --visualize

# 快速模式：单 partition 简单场景
beast2py quick \
    --alignment examples/primates.fasta \
    --tree-prior yule \
    --subst-model hky \
    --clock-model strict \
    --calibration-yaml examples/calibrations.yaml \
    --chain-length 1000000 \
    --name quick_analysis \
    --output output_quick.xml

# 验证生成的 XML（结构检查；加 --beast2 走 BEAST2 闸门）
beast2py validate --xml examples/output_basic.xml

# 生成 LaTeX 方法学描述
beast2py methods --config examples/config_basic.yaml --output methods.tex

# 生成分析指纹
beast2py fingerprint --config examples/config_basic.yaml --output fingerprint.json

# 生成 Snakemake/Nextflow 管道
beast2py pipeline --config examples/config_basic.yaml --type snakemake --output-dir pipeline

# 列出所有支持的模型
beast2py list-models
```

`generate` 与 `quick` 在没有 `--force` 时拒绝覆盖已存在的 `--output`；只要运行过的
验证闸门中有任何一个失败就返回非零，并且**不写出任何文件**（见[验证](#验证)）。

### API 使用

```python
from beast2py.api import Beast2Py

# 创建 API 实例
b2p = Beast2Py()

# 从 YAML 配置文件生成 XML。给出 output 时，XML 先走与 CLI 相同的验证闸门，
# 全部通过才写盘；目标文件已存在时需要 force=True。
xml_str = b2p.generate_xml("examples/config_basic.yaml", output="output_basic.xml")

# 运行校准诊断
report = b2p.diagnose("examples/config_basic.yaml", output_dir="diagnostics/", report_path="report.html")
for c in report.conflicts:
    print(f"[{c.severity}] {c.description}")

# 验证生成的 XML
result = b2p.validate_xml("output_basic.xml")
print(f"Valid: {result.is_valid}")

# 生成 LaTeX 方法学描述
methods = b2p.generate_methods("examples/config_basic.yaml", output="methods.tex")

# 生成分析指纹
fingerprint = b2p.generate_fingerprint("examples/config_basic.yaml")
print(f"Analysis fingerprint: {fingerprint}")

# 列出所有支持的模型
models = b2p.list_models()
print(models["substitution_models"])
```

`generate_xml()`、`quick_generate()` 与一体化 `generate()` 都接受 `force=False`。闸门失败时抛出
`ValueError`，并保持输出文件原样。

## 文档

详细的文档位于 `docs/` 目录中：

- **[使用手册](docs/user_manual_zh.md)** — 完整的安装、CLI/API、配置参考、诊断和故障排除手册
- **[架构指南](docs/architecture_zh.md)** — 分层架构和模块设计的详细描述
- **[教程](docs/tutorial_zh.md)** — 带有示例配置的分步教程
- **[校准指南](docs/calibration_guide_zh.md)** — 指定校准先验和来源元数据的指南
- **[XML 格式指南](docs/xml_format_zh.md)** — 生成的 BEAST2 XML 格式描述

每份手册还各有一节**严格校验**，逐条列出解析器拒绝的输入以及对应的报错原文。

## 示例配置

`examples/` 目录中共有 **19** 个完整配置。`tests/test_beast2_validation.py` 会为全部 19 个配置
生成 XML，并交给真实的 BEAST 2.7.8 `XMLParser` 检查（既解析，也调用 `initAndValidate` 初始化
模型）。它们实际覆盖的模型范围为：

- **替换模型：** 13 个支持模型中 9 个有示例——7 个核苷酸模型全部覆盖（JC69、HKY、TN93、GTR、SYM、TIM、TVM），氨基酸矩阵覆盖 WAG、JTT 两个；`dayhoff`、`blosum62`、`cprev`、`mtrev` 支持但无示例
- **时钟模型：** 全部 4 个（strict、UCLN、UCE、RLC）
- **树先验：** 8 个支持先验中 7 个有示例；`coalescent_constant` 无独立示例
- **校准分布：** 11 个支持分布中示例实际用到 3 个（lognormal、normal、uniform）

| 文件 | 场景 | 树先验 | 替换模型 | 时钟模型 | 校准数 | 特殊功能 |
|---|---|---|---|---|---|---|
| `config_basic.yaml` | 基础分析 | Yule | HKY+Γ4 | Strict | 1 | 显式 `parameter_priors` + 来源元数据；`output_basic.*` 的来源 |
| `config_advanced.yaml` | 高级多校准 | Birth-Death | GTR+Γ4+I / HKY | UCLN + 链接 | 3 | 密码子分区 + 位点/时钟链接 + 超先验 + `diagnostics` 块 + 茎群校准 |
| `config_subst_jc69.yaml` | 模型比较 JC69 | Yule | JC69+Γ4 | Strict | 1 | 每个候选模型单独一次运行 |
| `config_subst_sym.yaml` | 模型比较 SYM | Birth-Death | SYM+Γ4 | UCLN | 2 | `rates` 向量 + `frequencies.mode: uniform` |
| `config_subst_tim.yaml` | 模型比较 TIM | Yule | TIM+Γ4 | Strict | 1 | 具名速率参数 + `frequencies.mode: empirical` |
| `config_subst_models.yaml` | 模型比较 TVM | Yule | TVM+Γ4 | Strict | 1 | 五个具名速率参数 |
| `config_aa_model.yaml` | 氨基酸模型 | Birth-Death | WAG+Γ4 | UCLN | 2 | `aminoacid.fasta` + 带 offset 的 lognormal 校准 |
| `config_aa_jtt.yaml` | 氨基酸模型 | Yule | JTT+Γ4+I | Strict | 2 | 蛋白数据上的 Gamma+I |
| `config_tip_dates.yaml` | 末端日期 | Birth-Death | HKY | UCLN | 1 | `TraitSet`、`date-backward`、每个分类单元都有日期 |
| `config_calibrated_yule.yaml` | CalibratedYule | CalibratedYule | HKY+Γ4 | Strict | 3 | `calibration_method` 切换 + 茎群校准 |
| `config_multi_partition.yaml` | 多 partition | Birth-Death | HKY+Γ4 / GTR+Γ4 | UCLN + 链接 | 2 | `FilteredAlignment` 密码子分区、时钟链接、单一连锁拓扑 |
| `config_relaxed_clock.yaml` | 松弛对数正态时钟 | Birth-Death | GTR+Γ4 | UCLN | 2 | 带上下界的 `ucld_mean` / `ucld_stdev` |
| `config_relaxed_exponential.yaml` | 松弛指数时钟 | Yule | GTR+Γ4 | UCE | 2 | UCE 松弛时钟 |
| `config_rlc.yaml` | 随机局部时钟 | Yule | HKY+Γ4 | RLC | 1 | 指示变量与速率状态节点 |
| `config_skyline.yaml` | 贝叶斯天际线 | Bayesian Skyline | TN93+Γ4 | Strict | 1 | 天际线群体函数 |
| `config_coalescent_exponential.yaml` | 指数增长协方差 | Coalescent（指数） | HKY+Γ4 | Strict | 1 | 可正可负的 `growth_rate` |
| `config_ebsp.yaml` | EBSP | EBSP | HKY+Γ4 | Strict | 1 | 扩展贝叶斯天际线（BEAST 2.7 核心包内置） |
| `config_bd_skyline.yaml` | BD Skyline Serial | BD Skyline Serial | HKY+Γ4 | Strict | 1 | 需 [bdsky](https://github.com/BEAST2-Dev/bdsky) 插件包 |
| `config_nested_sampling.yaml` | Nested Sampling | Yule | HKY+Γ4 | Strict | 1 | NS MCMC（需 [NS](https://github.com/BEAST2-Dev/nested-sampling) 插件包） |

此外，`calibrations.yaml` 是给 `quick` 子命令使用的独立校准点文件。示例数据包括 `primates.fasta`（12 个分类单元 × 898 bp，由 BEAST2 自带示例比对 `Primates.nex` 转换而来）和 `aminoacid.fasta`（10 个分类单元 × 234 aa）。

## 仓库结构

- `beast2py/` — Beast2Py Python 包源代码（CLI、API、配置解析、XML 生成、诊断、可重复性、验证等模块）
- `examples/` — 19 个示例配置文件、示例比对数据（FASTA）和示例输出制品
- `tests/` — 测试套件（共收集 347 个测试：326 个单元/语义/发布完整性测试 + 21 个 BEAST2 集成验证测试）
- `beast2py/tools/` — 随包发布的无头 BEAST2 验证助手（`Beast2Validator.java` 源码，含预编译的 `classes/` 与 `launcher.jar`），由 `pyproject.toml` 以 `package-data` 随 wheel 打包发布；仓库根目录仅保留同名的薄包装 `beast2_validate.sh`
- `docs/` — 使用手册、架构指南、教程、校准指南和 XML 格式指南（中英文）
- `scripts/check_figure_export.py` — 三张 draw.io 图的导出后检查（它们没有可重绘的 Python 源）：读回已交付的 PDF 与 PNG，若导出件早于其 `.drawio` 源、字号低于 7.92 pt、字体超出 Helvetica 族、线宽低于 0.71 pt、或位图在页面自身宽度下不足 600 ppi，即报错退出（需 `pip install -e .[figures]`）
- `scripts/check_submission_placeholders.py` — 上传前检查：列出稿件或投稿信中仍遗留的全部未替换存档标识符（Zenodo DOI、Dryad DOI 与审稿链接、TreeBASE 登录号），并给出其段落或页码；只要还有未替换项即以非零码退出（PDF 检查需 `pip install -e .[figures]`）
- `beast2_validate.sh` — 无头 BEAST2 验证启动脚本（需要 `BEAST.base.jar` 与 JDK 17+；自动探测 `~/.beast/2.7/BEAST.base/*/lib/BEAST.base.jar`，也可设 `BEAST2_JAR`）
- `pyproject.toml` — 包构建与依赖配置

## 验证

`generate` 与 `quick` 依次运行三道闸门。只要运行过的闸门报告失败，命令就返回非零，并且
**不写出输出文件**。输出文件已存在时需要 `--force`，否则拒绝；同级 `.fingerprint.json` /
`.methods.tex` 将变为陈旧文件时，命令会给出警告。三道闸门是：

1. **结构检查** — XML 格式正确性与必需元素、id/idref 一致性、重复 ID、初值落在边界内、序列数据
   完整性；同时确认每个被估计参数都带 `Prior`、每个被估计状态节点都配有算子、没有算子挂在固定参数上
2. **校准冲突检测** — 嵌套支（含根校准）之间的时间一致性、分布重叠、单系性冲突
3. **BEAST2 检查（可选，`--beast2-validate`）** — 无头 BEAST2 `XMLParser` 在解析之后会对每个
   对象调用 `initAndValidate()`，因此检查的是模型能否**完成初始化**，不只是能否解析

**验证不能证明什么**：三道闸门都不能判断模型是否回答了你的生物学问题。BEAST2 能顺利初始化的
模型，仍然可能是错误的模型。

闸门 3 只在装有 BEAST2 与 JDK 17 的环境里才会运行。若你显式要求 `--beast2-validate`，而环境中
找不到 BEAST2，生成会以退出码 2 停止且不写任何文件——只读退出码的 CI 因此不会把“BEAST2 缺失”
当成通过。`--allow-unvalidated` 可以恢复“照样写出”的行为，此时不应把该 XML 称为“已通过 BEAST2
验证”。所有结论仅针对 BEAST2 **v2.7.8** 核实。

测试套件共收集 347 个测试（326 个单元/语义/发布完整性测试 + 21 个 BEAST2 集成测试）。集成测试对全部 19 个
示例配置、`quick` 路径以及仓库内的 `output_basic.xml` 生成 XML，并用真实的 BEAST 2.7.8 解析器
与模型初始化检查。当环境缺少 BEAST2 jar 或 JDK 17 时，这些测试直接跳过，不算失败。

## 可重复性

Beast2Py 生成以下可重复性制品：

- **分析指纹** — 确定性标识符 `B2P-{sha256(config)[:12]}-{version}`（例如 `B2P-ba18c26ed61f-0.1.0`）。其中**不含日期分量**，因此同一配置加同一比对数据，在任何时区、任何时刻重跑都得到相同的字节与相同的标识符。只改输出文件名不会改变指纹；序列内容、树设置、冠群/茎群标志、超先验、参数先验、末端日期、算子权重，任何一项变化都会改变指纹
- **数据摘要** — XML 注释里除标识符外还带一个独立的摘要：16 位十六进制，覆盖比对的**内容**。生成时刻的墙钟时间只出现在 `.fingerprint.json` 侧车文件里，绝不写进 XML
- **LaTeX 方法学描述** — 可直接用于论文的方法学段落
- **Snakemake/Nextflow 管道** — 端到端可重复工作流文件
- **校准来源** — 每个校准点的化石来源、参考文献 DOI 和校准类型

## 路线图

- 发布到 PyPI（`beast2py` 包名已核实可用）
- 支持 CladeAge 校准类型（基于 Birth-Death 模拟的校准密度，映射到 `cladeage.math.distributions.FossilPrior` XML）
- prior-only 日志的 ESS/HPD 自动解析与阈值告警（参考 Tracer 式 ACT/ESS 算法）
- BEAST2 运行编排（pybeast 式独立运行目录、种子管理、resume 快照）
- BEAST2 版本自动检查
- BEAST2 v2.6.x 命名空间兼容
- pirouette 式端到端实证测试（已知真树 → 模拟序列 → 注入冲突校准 → 检验冲突检测器灵敏度）
- 关注 BEAST 3 / LPhy 生态演进并适时适配

## 第三方许可

`beast2py/tools/classes/beast/pkgmgmt/` 下随包分发的编译类来自
[BEAST 2](https://github.com/CompEvol/beast2) v2.7.8，以 **GNU LGPL v2.1** 许可；
本包其余部分为 MIT。许可正文见 `beast2py/tools/LICENSE.BEAST2-LGPL-2.1.txt`，
来源、版本、是否修改以及获取完整对应源码的方式见 `beast2py/tools/NOTICE`。
`beast2py/tools/Beast2Validator.java` 随包分发其源码，以满足 LGPL-2.1 对该编译单元的要求。

## 许可证

本包为 MIT；随包捆绑的 BEAST2 包管理类为 LGPL-2.1（见 `beast2py/tools/NOTICE`）。SPDX 标识为 `MIT AND LGPL-2.1-only`。详见 [LICENSE](LICENSE)。

## 相关工具

BEAST2 生态中与 Beast2Py 相邻的开源工具（便于选型对比；亦见 BEAST2 官方博客《10 ways to generate BEAST XML》）：

| 工具 | 语言 | 定位 |
|---|---|---|
| BEAUti2（BEAST2 自带） | Java/JavaFX | 图形界面 XML 编辑器；v2.7+ 内置交互式方法学段落查看器与先验分布预览 |
| [beautier](https://github.com/ropensci/beautier) / [babette](https://github.com/ropensci/babette) | R | 程序化 XML 生成；babette 套件（beastier 运行、mauricer 插件管理、tracerer ESS/HPD 诊断）提供 R 一站式工作流 |
| [BEASTling](https://github.com/lmaurits/BEASTling) | Python | INI 配置 → BEAST2 XML（面向语言学/CLDF 数据；含校准与单系约束冲突检查） |
| [beast2-xml (acorg)](https://github.com/acorg/beast2-xml) | Python | 基于 BEAUti 模板的命令行 XML 生成器（FASTA/FASTQ 输入） |
| [BEASTmasteR](https://github.com/nmatzke/BEASTmasteR) | R | NEXUS + Excel 工作簿 → XML，面向化石 tip-dating / FBD |
| [LPhy](https://github.com/LinguaPhylo/linguaPhylo) / [LPhyBeast](https://github.com/LinguaPhylo/LPhyBeast) | Java | 概率模型规约语言 → BEAST XML（已面向 BEAST 3） |
| [pybeast](https://github.com/Wytamma/pybeast) / [beastiary](https://github.com/Wytamma/beastiary) | Python | BEAST2 运行封装（独立运行目录/种子/resume）与实时 MCMC trace 监控（ESS/HPD） |

Beast2Py 的差异点在于把以下能力整合为一条配置阶段流水线：YAML 声明式配置、三重校准冲突检测（时间一致性、分布重叠、单系性）、留一法先验敏感性 XML 的自动编排，以及校准来源（DOI）、分析指纹、方法学描述和管道生成。

## 引用

如果你在研究中使用 Beast2Py，请引用：

> Beast2Py: A Python framework for reproducible divergence time estimation with automated calibration prior specification, validation, and diagnostics.（手稿准备中）

## 作者

作者信息将在正式发布时补充。

## 资助

资助信息将在正式发布时补充。
