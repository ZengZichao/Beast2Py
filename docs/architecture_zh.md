# 架构设计

## 概述

Beast2Py 采用分层架构设计，注重可扩展性和可维护性。系统分为五层：CLI 入口层、配置管理层、领域模块层、数据模型层和工具支持层。

## 模块职责

| 模块 | 文件 | 职责 |
|------|------|------|
| CLI 入口 | `main.py` | 命令行参数解析、三道验证闸门、流程调度 |
| Python API | `api.py` | `Beast2Py` 门面类——与 CLI 完全同一条流水线，可从 Python 调用 |
| 配置管理 | `config.py` | YAML 配置解析、严格的键/布尔/值域校验、默认值填充 |
| 数据模型 | `models.py` | 所有数据类定义 |
| 序列读取 | `sequence.py` | 通过 biopython 读取 FASTA/NEXUS 文件 |
| Partition 管理 | `partition.py` | 多 partition 模型 link/unlink 管理 |
| 校准组装 | `calibration.py` | MRCAPrior、TaxonSet、分布 XML 构建 |
| 模型组装 | `model.py` | 树先验、替换模型、时钟模型、操作器、状态、初始化 XML 构建 |
| XML 生成 | `xml_writer.py` | 使用 ElementTree 组装完整 XML，并维护操作器/状态节点不变式 |
| 校准诊断 | `diagnostics.py` | 冲突检测、敏感性分析、可视化（核心创新） |
| 可重复性 | `reproducibility.py` | 分析指纹、方法学描述、管道生成（核心创新） |
| 验证器 | `validator.py` | XML 结构与语义校验、id/idref 一致性、BEAST2 闸门 |
| 工具函数 | `utils.py` | ID 生成、XML 序列化、分布计算 |
| 模型注册表 | `registry.py` | 模型名称到 XML spec 的映射、可扩展注册 |
| BEAST2 校验启动脚本 | `beast2_validate.sh` | 无头 BEAST2 验证入口，发布在**包内**（仓库根只保留同名的薄包装） |
| Java 校验助手 | `tools/` | `Beast2Validator.java` 源文件、预编译的 `classes/` 与 `launcher.jar`，作为 package-data 发布；`NOTICE` 与 `LICENSE.BEAST2-LGPL-2.1.txt` 记录这些捆绑类的来源与所承载的 LGPL-2.1 条款 |

## 数据流

```
1. 用户输入
   ├── 序列文件（FASTA/NEXUS，支持多个）
   └── YAML 配置文件（含校准信息 + provenance）
        │
        ▼
2. 配置解析 (config.py) —— 在任何东西被构建之前，严格、失败即关闭
   ├── 序列读取 (sequence.py) → List[Alignment]
   ├── Partition 管理 (partition.py) → List[Partition]（alignment + 模型配置 + link/unlink）
   └── 在此被拒绝：未知的键、含混的布尔值、越界的先验、非法的 FilteredAlignment
       表达式、`tree: separate`，以及被*硬*边界违背的 provenance 年龄记录（*软*边界只告警）
   → BEASTConfig 对象（partition 列表 + 校准点列表 + MCMC + diagnostics + provenance）
        │
        ▼
3. 模型组装 (model.py)
   → 各 partition 的树先验 XML + 替换模型 XML + 时钟模型 XML，共享同一棵连锁拓扑；
     每个 `estimate: false` 的参数都被*内联*写出，而不进入状态空间
        │
        ▼
4. 校准组装 (calibration.py)
   → MRCAPrior XML 元素 + TaxonSet 元素（末端日期还有 TraitSet）
        │
        ▼
5. XML 生成 (xml_writer.py) —— 仅在内存中
   → 使用 xml.etree.ElementTree 构建完整 XML；写入器同时断言每个被估计状态节点
     都有算子、且没有算子挂在被固定的参数上
        │
        ▼
6. 验证闸门 (main.py::_validate_and_write) —— 所有运行过的闸门都通过之前不写任何文件；
   任一失败都以非零码退出
   ├── 闸门 1 结构与语义校验 (validator.py)
   ├── 闸门 2 校准冲突检测 (diagnostics.py) ← 核心创新
   │     时间一致性、分布重叠、单系性冲突
   └── 闸门 3 可选的 BEAST2 解析 + initAndValidate (beast2_validate.sh + tools/)
        │
        ▼
7. 输出（仅在闸门通过后）
   ├── BEAST2 XML 文件
   └── 可重复性打包 (reproducibility.py) ← 核心创新
       ├── 分析指纹（配置哈希 + 数据哈希 + 工具版本，不含日期）
       ├── 方法学描述（LaTeX 片段）
       └── Snakemake/Nextflow 管道文件（独立的 `pipeline` 命令）

用同一批构件实现的独立命令：
   diagnose (diagnostics.py) → 冲突报告（HTML/文本）、N+1 个仅先验采样 XML、
   校准分布密度图
   methods / fingerprint / pipeline → 单独产出各类可重复性制品
```

## 模型注册表模式

`registry.py` 中的 `ModelRegistry` 类实现了注册表模式，把配置中的模型名称映射到对应的 BEAST2 XML spec 值。这带来三个好处：

- **可扩展性**：用户可注册自定义模型，无需修改核心代码
- **一致性**：所有模型 spec 集中定义
- **验证**：未知模型名称在解析时即可检测

## 关键设计决策

1. **使用 ElementTree 生成 XML**：使用 `xml.etree.ElementTree` 构建 XML DOM 树后序列化，保证 XML 转义和格式一致性

2. **YAML 配置**：声明式 YAML 配置文件，便于可重复性和版本控制

3. **失败即关闭的校验**：解析器把每个配置段都对照键白名单校验，拒绝含混的布尔值、越界的先验以及生成器无法兑现的输入，而不是静默地填入默认值。三道验证闸门在**写出任何文件之前**运行，任一失败都以非零码退出。`--force` 只放开其中唯一安全的一项：替换一个已存在的输出文件

4. **被估计参数与被固定参数**：生成器把标注 `estimate: false` 的参数写成使用者所在元素的内联 `<parameter estimate="false">` 子元素，因此该参数既不进入 `<state>`，也不获得算子，也不产生轨迹列（`model.py::_state_node`）。写入器随后双向守护这一不变式：被估计的状态节点缺少算子、算子挂在被固定的参数上，这两种情况都会报错

5. **检查点间隔随链长伸缩**：`mcmc.store_every` 缺省为 `max(1, chain_length // 10000)`（`config.py`），因此任何运行都能得到大约十个检查点，而不是只有最末一个。生成器刻意**不**把随机数 `seed` 写进 XML：BEAST2 的 `MCMC` 没有 `seed` 输入，种子只能通过命令行参数 `beast -seed N` 传入

6. **Partition ID 命名空间**：所有参数 ID 使用 partition ID 前缀（如 `gene1.hky.kappa`），支持多 partition

7. **单一连锁拓扑**：所有 partition 都装配在同一棵共享树上。解析阶段直接拒绝 `tree: separate` 与逐分区的 tree id：若静默退化为连锁基因树，本意为 *BEAST 风格的分析就变成了完全连锁的分析

8. **尺度无关的诊断**：分布重叠判据默认是无量纲的——用 1-Wasserstein 距离除以两个先验 95% 区间的平均宽度，阈值为 0.15。旧有的绝对值判据（2.0 个时间单位）仍可通过 `diagnostics.overlap_measure` 选用。该检查会跳过没有定义统计量的分布（OneOnX）

9. **自包含的 BEAST2 闸门**：无头验证器位于包内（`beast2py/beast2_validate.sh`、`beast2py/tools/`），因此 `pip install .` 就自带可用的闸门 3。脚本按一批候选路径搜索 JDK 17 或更新版本（`BEAST2_JAVA_CANDIDATES`、`JAVA_HOME`、`/usr/libexec/java_home`、Homebrew 与系统 JVM 目录，最后才是 `PATH` 上的 `java`），并支持 `BEAST2_JAR`、`BEAST2_SKIP_JAVA_HOME_TOOL`、`BEAST2_SKIP_JDK_PROBES` 与 `BEAST2_VALIDATE_TRACE`

10. **不含日期的指纹**：分析指纹为 `B2P-{hash12}-{version} | Data: {hash16}`，即配置摘要加上比对序列的摘要。同一配置与同一数据在任何机器、任何日期都得到同一标识符；墙钟生成时间只存在于 JSON 侧车文件

11. **BEAST2 v2.7.x 命名空间**：使用 `beast.base.*` 包名，已针对 BEAST 2.7.8 核实——该版本的 class 文件在运行时需要 JDK 17 或更新版本
