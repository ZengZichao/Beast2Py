# 校准先验最佳实践指南

## 引言

节点校准是分歧时间估计的关键步骤。不当的校准设定会让后验估计产生偏差，甚至出现“先验支配”。本指南给出在 Beast2Py 中指定校准先验的最佳实践。

## 关键原则

### 1. 使用适当的分布

- **正态分布**：适用于不确定性对称的精确校准
- **对数正态分布**：适用于右偏不确定性（化石数据中常见）
- **均匀分布**：仅适用于硬边界约束
- **指数分布**：适用于只有最小年龄约束的情况

Beast2Py 共接受 11 种分布类型。上述四种用来描述节点年龄，在生物学上通常说得通。其余为 `gamma`、`beta`、`laplace`、`inverse_gamma`、`poisson`、`chi_square`（可参数化的形状），以及 `one_on_x`（`f(x) = 1/x`）。后者是**非正常密度**（improper）：既没有均值，也没有 95% 区间，因此 Beast2Py 把它的统计量报为“未定义”，并把它排除在重叠检查之外。只有当参数的范围由别的东西界定时，才适合用 `one_on_x` 作为刻意平坦的先验。

### 2. 软边界与硬边界

- **软边界**（推荐）：允许后验超出先验范围，让数据参与估计
- **硬边界**：严格约束后验；仅在化石证据确定时使用

每个校准点用 `provenance.calibration_type` 记录这一选择（自由字符串；除 `"hard"` 之外的任何值都按软边界处理，省略时默认为**软边界**）。Beast2Py 还会*兑现*你记录的声明：解析阶段把先验的 95% 区间与 `provenance.original_age_min` / `original_age_max` 比较，于是

- 校准声明为**硬边界**，而先验整体落在所记录的化石范围之外，即为**错误**——流程在写出任何文件之前终止；
- 同样情形下，**软边界**的校准只给出**警告**，因为软边界本就把大部分质量放在化石最小值之上（让 lognormal 的 2.5% 分位数落在硬最小值处，正是 BEAST2 的惯例）。

### 3. 记录来源

始终为每个校准点指定 `provenance` 块：

```yaml
calibrations:
  - name: "homininaMRCA"
    taxa: ["Homo_sapiens", "Pan"]
    distribution:
      type: "lognormal"
      parameters: {M: 0.35, S: 0.4}
      offset: 6.5
    provenance:
      source: "fossil"           # fossil / biogeographic / secondary / other
      reference: "DOI:10.1038/nature12323"
      calibration_type: "soft"    # soft / hard
      original_age_min: 7.0
      original_age_max: 9.5
      notes: "Sahelanthropus tchadensis, Toros-Menalla 7.24 Ma; soft minimum"
```

`provenance` 恰好接受 `source`、`reference`、`calibration_type`、`original_age_min`、`original_age_max` 与 `notes` 六个键，解析阶段会拒绝其他任何键；`original_age_min` 不得大于 `original_age_max`。

### 4. 检查冲突

使用 `diagnose` 命令检测：

- **时间不一致**：嵌套 clade 的父节点比子节点更年轻。父 clade 的 HPD 上界比子 clade 的 HPD 下界更年轻时判为错误，只有均值冲突时判为警告。根校准（`taxa: null`）充当其余每个 clade 的父节点参与比较
- **分布重叠**：两个非嵌套校准的先验几乎可以互换。默认判据是**无量纲的**：两个先验之间的 1-Wasserstein 距离，除以它们 95% 区间的平均宽度，低于 `0.15` 时告警。固定距离会随节点深度增长，因此在浅节点上容易触发，在深节点上几乎不会触发。旧有的绝对判据（Wasserstein 距离低于 `2.0` 个时间单位）仍然可用：把同一 `diagnostics` 段里的 `overlap_measure` 设为 `absolute`（或 `both`），并用 `overlap_threshold` / `relative_overlap_threshold` 调整切点。`OneOnX` 没有定义的均值与区间，重叠检查会跳过它
- **单系性冲突**：有交集但不嵌套的分类子集

同一套冲突检测会作为 `generate` 与 `quick` 的闸门 2 自动运行：出现 `error` 级别的冲突时，闸门不写出任何文件，命令以非零码退出。

### 5. 运行敏感性分析

始终评估每个校准点的边际影响：

```bash
beast2py diagnose --config config.yaml --report report.html --sensitivity
```

此命令生成 N+1 个仅先验采样 XML（N 为校准点数量）：一个包含完整校准集，外加 N 个留一法版本。它们落在 `<报告目录>/diagnostics/sensitivity/` 下，或落在 `--sensitivity-xml-dir` 所指目录的 `sensitivity/` 子目录里。用 BEAST2 运行每个 XML，然后比较树高的后验。短链即可：1,000,000 代，即 `mcmc.chain_length: 1000000`。

## 常见陷阱

### 1. 先验支配

**问题**：校准先验压倒序列数据，使后验与先验无法区分。

**解决方案**：把仅先验采样的结果与后验作对比。如果两者相同，说明先验的信息量过大。模型过度受约束时，不额外运行一次也能看出来：轨迹只是跟着先验的 HPD 走。

### 2. 冲突校准

**问题**：两个校准为嵌套 clade 指定了不兼容的时间范围。

**解决方案**：运行 `diagnose` 检测时间冲突，然后调整校准，让父 clade 的范围比子 clade 更老。注意：闸门 2 报出 `error` 级别的冲突时，`generate` 不会写出任何文件，因此退出码为 0 的流水线已经通过了这一检查。

### 3. 过窄先验

**问题**：使用很窄的正态分布（`sigma` 很小）会人为地约束后验。

**解决方案**：对化石校准使用更宽的分布或对数正态分布，因为化石校准通常具有不对称的不确定性。

### 4. 来源记录与实际先验不一致

**问题**：记录的化石年龄范围与实际使用的先验是两个不同的年龄。BEAST2 只看得到分布，看不见这类记录错误。

**解决方案**：如实填写 `original_age_min` / `original_age_max`。Beast2Py 会在解析阶段把这两个范围与先验的 95% 区间比较：硬边界不一致判为错误，软边界不一致判为警告，提醒你在方法学里把两者一起交代清楚。

## 参考文献

- Heled & Drummond (2012) Syst Biol 61: 138～149, DOI:10.1093/sysbio/syr087 — CalibratedYule 和 MRCA 校准方法
- Rieux & Balloux (2016) Mol Ecol 25: 1911～1924, DOI:10.1111/mec.13586 — 校准最佳实践
- Ho & Phillips (2009) Syst Biol 58: 367～380, DOI:10.1093/sysbio/syp035 — 分子钟校准
- Parham et al. (2012) Syst Biol 61: 346～359, DOI:10.1093/sysbio/syr107 — 化石校准最佳实践
