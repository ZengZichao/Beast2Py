# BEAST2 XML 格式指南

[English](xml_format_en.md) | 中文

## 概述

Beast2Py 生成兼容 BEAST2 v2.7.x 的 XML 文件，使用 `xml.etree.ElementTree` 构建，并已在 BEAST 2.7.7 运行时上核实（其中 `BEAST.base` 包升级到 2.7.8）。本文档描述工具生成的 XML 结构。工具在内存中组装 XML，验证闸门（结构与语义检查、校准冲突检测，以及按需开启的 BEAST2 解析 + `initAndValidate()`）全部通过后，才把结果写到磁盘。因此磁盘上的文件一定具有下述结构。

## XML 结构

```xml
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<beast version="2.0" namespace="...">
    <!-- 1. 序列比对数据 -->
    <data id="alignment" dataType="nucleotide">
        <sequence taxon="human">ATGGCAAAT...</sequence>
        ...
    </data>

    <!-- 可选：密码子分区的 FilteredAlignment -->
    <data data="@alignment" dataType="nucleotide" filter="3::3" id="thirdSites"
          spec="beast.base.evolution.alignment.FilteredAlignment"/>

    <!-- 2. 校准用的分类单元集合 -->
    <!-- 实现约定：首次出现的分类单元用 id= 定义，之后引用用 idref=。
         全部用 idref，当不存在顶层的分类单元定义时会导致 BEAST2 报
         "Could not find object associated with idref ..." -->
    <taxonset spec="TaxonSet" id="humanChimpMRCA.taxonset">
        <taxon spec="Taxon" id="human"/><taxon spec="Taxon" idref="chimp"/>
    </taxonset>

    <!-- 3. 后验分布 -->
    <distribution id="posterior" spec="CompoundDistribution">
        <distribution id="prior" spec="CompoundDistribution">
            <!-- 树先验 -->
            <!-- MRCA 校准先验 -->
            <!-- 参数先验 -->
        </distribution>
        <distribution id="likelihood" spec="CompoundDistribution">
            <!-- 每个分区的 TreeLikelihood -->
        </distribution>
    </distribution>

    <!-- 4. MCMC 运行配置 -->
    <run id="mcmc" spec="MCMC" chainLength="10000000" preBurnin="0" storeEvery="1000">
        <init spec="beast.base.evolution.tree.coalescent.RandomTree" id="randomTree"
              initial="@tree" taxa="@alignment" rootHeight="10.275">...</init>
        <state id="state" storeEvery="1000">...</state>
        <distribution idref="posterior"/>
        <operator>...</operator>
        <logger id="tracelog" .../>
        <logger id="treelog" .../>
        <logger id="screenlog" .../>
    </run>
</beast>
```

### `<state>` 元素只容纳被估计的参数

用户标注 `estimate: false` 的参数永远不进入 `<state>`。生成器把这类参数*内联*写成引用它的元素的子元素，因此它不获得算子，也不产生轨迹列：

```xml
<!-- 被估计：一个状态节点，由替换模型的 kappa 属性引用 -->
<parameter id="alignment.hky.kappa" name="stateNode" value="2" dimension="1"
          estimate="true" lower="0"/>
<substModel spec="HKY" id="alignment.hky" kappa="@alignment.hky.kappa">...</substModel>

<!-- 被固定：内联参数，无状态节点、无算子、不被记录 -->
<substModel spec="HKY" id="alignment.hky">
    <parameter id="alignment.hky.kappa" name="kappa" value="2" estimate="false"
               dimension="1" lower="0"/>
</substModel>
```

写入器双向守护这一不变式，宁可报错，也不产出静默出错的模型：每个被估计的状态节点都必须配有算子（否则 BEAST2 只是给个警告，然后把该参数钉在初值上），算子也只能挂在被估计的参数上。

### `<run>` 元素的属性

`chainLength`、`preBurnin`、`storeEvery` 写在 `<run>` 上（`storeEvery` 同时也写在 `<state>` 上）。`storeEvery` 计量的是**已记录的样本数**，不是代数；除非 `mcmc.store_every` 另有指定，它默认为 `max(1, chain_length // 10000)`。XML 中刻意**没有 `seed` 属性**：BEAST2 的 `MCMC` 没有 `seed` 输入，写出一个会让 XML 无法解析——随机数种子是命令行参数 `beast -seed N`。日志文件名里仍可以出现 `$(seed)` 宏，由 BEAST2 在运行时展开。

当 `mcmc.type: nested_sampling` 时，`<run>` 元素改用附加包自己的引擎：

```xml
<run id="mcmc" spec="nestedsampling.gss.NS" chainLength="10000000" preBurnin="0"
     storeEvery="10000" particleCount="1" subChainLength="10000">
```

## 命名空间

对于 BEAST2 v2.7.x，命名空间使用 `beast.base.*` 包名：

```
beast.base.evolution.alignment:beast.base.evolution.speciation:beast.pkgmgmt:beast.base.core:...
```

## 多分区支持

生成器为每个分区写出一个独立的 `<data>` 元素，以及一个独立的似然节点 `<distribution id="gene1.treeLikelihood" spec="TreeLikelihood">`：其 `id` 模式是 `{partition}.treeLikelihood`（`gene1`/`gene2` 只是示意的分区名；随附的多分区示例用的是密码子分区 `codon12` 与 `codon3`）。并不存在 `<TreeLikelihood>` 这个元素标签，该类名只出现在 `spec` 属性里，而且写出器对每个属性都用双引号：

```xml
<data id="gene1" dataType="nucleotide">...</data>
<data id="gene2" dataType="nucleotide">...</data>

<distribution id="likelihood" spec="CompoundDistribution">
    <distribution id="gene1.treeLikelihood" spec="TreeLikelihood" data="@gene1" tree="@tree">
        <siteModel id="gene1.siteModel" spec="SiteModel">...</siteModel>
        <branchRateModel id="gene1.clockModel"
                         spec="beast.base.evolution.branchratemodel.UCRelaxedClockModel"
                         normalize="true">...</branchRateModel>
    </distribution>
    <distribution id="gene2.treeLikelihood" spec="TreeLikelihood" data="@gene2" tree="@tree">
        <siteModel idref="gene1.siteModel"/>           <!-- 共享位点模型 -->
        <branchRateModel idref="gene1.clockModel"/>    <!-- 共享时钟 -->
    </distribution>
</distribution>
```

共享的时钟只写出一次，并带上自己的 `id`，其余分区用 `idref` 指向它（共享位点模型同理）。所有分区都装配在**同一棵** `tree="@tree"` 上：逐分区拓扑尚未实现，解析器会直接拒绝 `tree: separate`，而不是静默地把各分区的树链接起来，因此产出文件里绝不会出现多个 `<tree>` 状态节点。

## 校准分布

所有分布支持 `offset` 属性（默认 0.0）。下表 `spec` 列给出类名；实际写入 XML 的是全限定形式，例如 `beast.base.inference.distribution.LogNormalDistributionModel`。

| 分布 | spec | 参数 |
|------|------|------|
| Normal | `Normal` | mean、sigma |
| LogNormal | `LogNormalDistributionModel` | M、S、meanInRealSpace |
| Uniform | `Uniform` | lower、upper |
| Exponential | `Exponential` | mean |
| Gamma | `Gamma` | alpha、beta |
| Beta | `Beta` | alpha、beta |
| OneOnX | `OneOnX` | （无） |
| Laplace | `LaplaceDistribution` | mu、scale |
| InverseGamma | `InverseGamma` | alpha、beta |
| Poisson | `Poisson` | lambda |
| ChiSquare | `ChiSquare` | df（配置键为 `dof`） |

`OneOnX`（`f(x) = 1/x`）是非正常密度：均值与 95% 区间都没有定义，因此诊断模块把它的统计量报为“未定义”，并把涉及它的配对排除在时间顺序与区间比较之外，将该对报告为 `unchecked`。

## 树先验

| 树先验 | spec |
|------|------|
| Yule | `YuleModel` |
| CalibratedYule | `CalibratedYuleModel` |
| Birth-Death | `BirthDeathGernhard08Model` |
| 溯祖（恒定种群大小） | `Coalescent` + `ConstantPopulation` |
| 溯祖（指数增长） | `Coalescent` + `ExponentialGrowth` |
| Bayesian Skyline | `BayesianSkyline` |
| EBSP | `Coalescent` + `CompoundPopulationFunction`（BEAST 2.7 核心包内置） |
| BD Skyline Serial | `bdsky.evolution.speciation.BirthDeathSkylineModel`（需要 bdsky 附加包） |

## 时钟模型

| 时钟模型 | spec |
|------|------|
| 严格时钟 | `StrictClockModel` |
| UC 松弛对数正态 | `UCRelaxedClockModel` |
| UC 松弛指数 | `UCRelaxedClockModel` |
| 随机局部时钟 | `RandomLocalClockModel` |

UCLN 与 UCE 共用 `UCRelaxedClockModel` 这个 spec，区别只在分支速率分布：前者是 `M = 1`、`S = ucld.stdev` 的 `LogNormalDistributionModel`，后者是一个 `Exponential`，它的 `mean` 是写死的字面量 `<parameter name="mean" value="1.0"/>`；两种情况下，用户侧的 `ucld.mean` 都不接进这个分布，而是作为 `clock.rate` 输入挂在 `UCRelaxedClockModel` 元素上。只有 UCLN 与 UCE 这两个松弛型的 `UCRelaxedClockModel` 元素带 `normalize="true"`，与 BEAUti 的参考版式一致，从而保持速率与节点年龄可识别；随机局部时钟产出 `RandomLocalClockModel`，带 `ratesAreMultipliers="false"` 且**没有** `normalize` 属性，它还会额外引入两个专属状态节点 `rates` 与 `indicators`。

## 分析指纹

工具在 XML 注释中嵌入分析指纹：

```xml
<!-- Analysis Fingerprint: B2P-3f43bfab58c3-0.2.0 | Data: 72a8314e0bde0640 -->
```

格式：`B2P-{config_hash_12chars}-{tool_version} | Data: {data_hash_16chars}`。

`config_hash` 是序列化配置的 SHA-256 摘要前 12 位十六进制字符；`data_hash` 则把每个比对自身的 id 与内容摘要（分类单元名加序列字符，而非文件路径）折进一个 SHA-256 摘要，取其前 16 位十六进制字符。

该标识符刻意不含日历日期、时区，也不含输出文件名，因此同一配置与同一比对在任何机器、任何日期都得到同一指纹。改动任何模型组件、校准、MCMC 设置或序列数据，指纹都会随之改变；唯一有文档说明的例外是 `mcmc.threads`，它只调度运行、不改变后验，因此被排除在摘要之外。

生成时间只单独记录在 JSON 侧车文件里，绝不进入标识符。
