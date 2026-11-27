# 一作手册 · 从组织切片到可检验的结论

## 现场先看什么

打开 index.html，选一个实测公开队列。先读供者数和分组，再看采集位置，最后才看基因或离子的显著比例。灰点表示缺少审计结果。合成夹具可从来源菜单切换，每个夹具都提供构造过程和种子。网页本地打开即可运行，滑块对应的解析模型与仓库的队列模拟分别有说明。

## 看懂五个模块

**输入契约（schema.py）**　矩阵的一行可以是一个 spot、一个细胞或一个像素。每行的坐标须已校准到微米。每张切片归属于一个生物学个体，同一人多张切片共享 patient_id。若病人与对照共用载玻片，要留下真实 slide_id。同一患者出现在两个条件时，独立两组检验的交换单位不成立，当前程序会要求另行建立配对设计。

**提问一：实验能分辨多小的 p 值？**　design_audit.exact_floor 枚举思想是固定两组人数，在患者之间重新发标签。每个合法分法占同样大的概率格子。本工具把两个单侧尾部中较小的概率乘二，极端观测也受这个格子宽度限制。BH 看排序后的一整列 p 值，许多同样小的值能一起过线。最低可行秩描述分辨率；检验功效还取决于效应和方差。

**提问二：疾病标签和仪器设置能分开吗？**　design_audit.design_audit 把患者在多个片号/批次中的总权重限制为一，用归一化互信息衡量标签携带的信息。小型独立列联表还报告 Fisher 检验。线性规划检查两组的二维位置能否被直线完全分开。仅有切片本地坐标时，边界框中心只能提示构图或切割差异，需要结合共享载物台参考系解释。

**提问三：细胞很多，证据增加多少？**　pseudorep.pseudorep_audit 对每个特征先做合并细胞的秩和检验，然后把每张切片汇总、对同一患者等权平均，用患者层面的组标签置换。patient_test 在每个分法下重新估计跨特征的方差先验，将不稳定的个别方差向共同尺度收缩。公开队列的真值未知，两个 BH 调用率之差可以度量分析单位敏感性，不能直接识别每一个假发现。

**提问四：同一物理尺度的结构能被组成解释多少？**　spectral_audit.build_grid 将点汇入占据的方形控制体，使用共享面的通量构造扩散算子。频带是两个扩散核的差，经同一切比雪夫递推施加给 X 和 C。设某一频带的响应为 Y，取滤波后的 C 的数值列空间，得到正交基 Q。投影能量占比为 ‖QᵀY‖² / ‖Y‖²。零能量时该比值没有定义。解释占比是一项关联量，成分也可能处在真正的生物机制链上。

**提问五：白噪声能造成多大的高频变化？**　spectral_audit 使用相邻与第二距离的半方差估计短程截距，并按经过相同滤波与投影的噪声传递量修正细尺度能量。两段相邻边和跨越两段的远边来自同一直线路径。匹配路径截距等于二阶差分平方均值除以六，因而非负。患者级检验保留有符号的修正能量，包括估计噪声大于原始能量时的负值。对数归一化输入只支持变换尺度上的描述。短程突变、分割边缘与组织断裂也会贡献近邻差异，因此将触发项称为 nugget 敏感特征。

## 二十个追问与回答

1. **一个患者有更多切片，会得到更大权重吗？** 患者均值检验先在每张切片计算均值，再在患者内等权汇总。若按细胞池化，采样更密或组织更大的患者会贡献更大。
2. **精确置换何时有统计意义？** 在零假设下，给定组大小后患者标签可交换。分组与批次完全重合时，重新贴标签仍可计算参考分布，但疾病的单独效应缺少识别依据。
3. **不平衡组能用同一个地板公式吗？** 本实现选用双倍较小尾概率，两组不平衡时仍由同一组合数给出理论下界。绝对统计量双侧尾检验可能有不同下界。
4. **为什么不用首位 BH 临界值判断全组零功效？** BH 是逐步升阶程序。若足够多的特征并列在下界，较后秩可以达到更宽松的门槛。
5. **Monte Carlo 到底有多精确？** 超过穷举预算时，使用固定种子、有限抽样并对观测排列加一。报告写出模式、抽样数和理论分法数，两种分辨率应分别解读。
6. **经验贝叶斯会不会使用原标签形成优势？** 每次置换都重新拟合方差先验，使统计量以同样方式对待每个分法。小面板或尺度悬殊的离子仍会影响先验质量。
7. **细胞秩和和患者均值检验在测同一个量吗？** 统计量不同，前者衡量池化分布的秩差，后者比较患者平均水平。分布形状、测序深度和加权也能造成落差，报告保留两套结果。
8. **载玻片一人一张为什么也会触发？** 将每张片设为固定效应时可重建每个患者标签。独立随机分片有助于跨片推断，但需要随机化记录，单靠归档片号无法验证它。
9. **AUC 很高能证明混杂？** 在有限样本中，训练内 AUC 偏乐观。它描述当前坐标支持的重叠程度，交叉验证及设计记录可用于补充验证。
10. **把像素当作微米会怎样？** 热核的物理尺度会系统偏移。输入必须校准；缩放坐标与波长后能量占比应保持一致，性质测试覆盖该变换。
11. **数据中有组织孔洞怎么办？** 空控制体不参与图，分离的连通分量独立扩散。格宽以下的微小空洞可能被离散化合并。
12. **奇异的 C 会破坏回归吗？** 在滤波后用数值秩提取列空间，冗余列不会提供额外方向。若列空间几乎饱和，应结合独立注释与样本外验证。
13. **没有细胞类型比例，绿灯是否合理？** 成分维度要显示未测或部分可测。皮层区带指示量支持解释解剖区带的贡献，其范围与细胞比例不同。
14. **MSI 的 TIC 是否引入自相关？** TIC 包含目标离子，强峰会影响自身协变量。高解释率应结合留一离子 TIC、内标或外部采集漂移指标重新评估。
15. **nugget 估计是否依赖模型？** 它依赖近距离结构可外推、局部噪声近似白噪声及滤波噪声传递的近似。真实细尺度组织边界可能被过度修正。
16. **如何看待 nnSVG 或 SPARK-X 比 DESeq2 调用更多？** 前两者的空间变异与后者的组间表达差异有不同零假设。比较表按估计目标列出方法，仅在共同特征和适用终点内比较。
17. **缺包和无调用怎样区别？** 成功运行且没有特征过线意味着调用率为零；不适用、缺依赖和执行失败保留状态以及空数值。运行版本、时间与 RSS 随每行发布。
18. **百万细胞测试能代表任意机器吗？** 应汇报矩阵形状、稀疏度、每患者分块、选定特征、线程和峰值 RSS。增加稠密特征或更细网格需单独测量。
19. **为什么相图没有统一的最少患者数？** 检出率取决于预设效应、离散度、批次和空间场。表格给出同时满足检出率下置信界和零效应拒绝率上置信界的最小已测试每组患者数，完全批次对齐的设计留空。相图中的原始端点是未投影的频带能量，调整后同时扣除成分与估计噪声。D5 报告的扣噪前端点已完成成分投影，两者不能直接做数值对照。
20. **论文宣称与归档数据不同，如何下结论？** 检查模态、排除标准、训练与验证队列、组织标本与独立患者。保留原文位置和范围，交由人工复核；一个子集的结果不能推论整篇论文。

## 三分钟展示 · 中文

**第一分钟，界定问题。** 打开队列总览，选真实的 Visium 审计。说明灰色格子的含义，点击论文和存档出处。指出患者、切片、细胞三层单位。用精确检验地板解释为什么即便每个患者多测很多位置，独立标签分法数仍由患者数决定。

**第二分钟，展示证据。** 比较池化与患者级 BH 调用率，说明两个检验的 estimand。查看三个波长的解释占比，并说明所用 C 的生物学含义。打开采集坐标，报告显示参考系；若有共线警告，说明哪些效应无法独立分辨。切换到固定合成陷阱，把已知设计同公开案例的未知真值联系起来。

**第三分钟，说明边界和复现。** 在沙盒中保持患者数固定，仅增加每位患者的细胞数，观察解析误报曲线。进入实测 Benchmark，指出方法名称、终点和资源成本。最后说明如何从输入清单、SHA-256、结果 JSON 追到图表，展示 claim-check 的原文证据定位，并邀请读者用同一工具审计自己的设计。

## Three-minute demonstration · English

**Minute one.** Open a measured cohort and identify its archive, biological units and sections. Explain the two-sided permutation floor. Extra positions help estimate an individual tissue profile, while independent donors determine how many label assignments are available. Gray entries identify catalog records awaiting a measured audit.

**Minute two.** Compare pooled and patient-level discovery rates. Describe how rank distributions and patient means weight tissues differently. Read the three physical wavelengths and explain which covariates enter the projection. Inspect the acquisition coordinate frame and the observed overlap. Switch to a synthetic fixture to show what can be checked when the data-generating mechanism is known.

**Minute three.** Hold donor count fixed and increase cells per donor in the sandbox. Contrast the analytic illustration with the measured simulation grid. Show the recorded estimand, runtime and memory of an applicable comparator. Finish by tracing one figure from a checksummed input through a strict JSON result to its automatically generated manuscript value. The interpretation concerns reproducibility, resolution and sensitivity under the stated design.

