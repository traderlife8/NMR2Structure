# RDKit 功能总结

> 本文档参考 [RDKit 完整函数手册](doc/RDKit_manual.md) 提炼，面向功能选型、API 定位和日常化学信息学工作流。
>
> 本机依据：`rdkit 2026.03.4`、Python 3.13.5。完整函数签名、枚举、底层方法和可运行示例请查阅 [doc/RDKit_manual.md](doc/RDKit_manual.md)。
>
> 本文档不重复代码示例；所有 API 细节通过链接指向完整手册。

## 一、库概览

RDKit 是基于 C++ 化学核心的开源化学信息学工具包，通过 Python 提供分子表示、结构解析、属性计算、结构检索、指纹相似性、反应处理、2D/3D 建模、药效团分析、分子绘图和数据处理能力。

RDKit 的核心抽象是分子图：原子是节点，化学键是边，分子可以附带立体信息、多个二维或三维构象、属性、SGroup 和查询规则。它既能表示确定结构，也能表示 SMARTS 查询模板和反应模板。

### 核心功能地图

| 功能域 | 能力 | 主要模块 |
|---|---|---|
| 分子对象 | Mol、Atom、Bond、Conformer、RingInfo、SGroup | `rdchem` |
| 结构格式 | SMILES、SMARTS、InChI、MOL、SDF、PDB、CML、CDXML | `rdmolfiles` |
| 结构清洗 | 价态、芳香性、环、氢、立体和标准化 | `rdmolops`、`MolStandardize` |
| 结构搜索 | SMARTS、子结构库、合成空间、R 基团 | `rdqueries`、`rdSubstructLibrary` |
| 描述符 | 质量、LogP、TPSA、氢键、环、拓扑、3D | `Descriptors`、`rdMolDescriptors` |
| 指纹 | Morgan、RDKit、MACCS、AtomPair、Torsion、Avalon、MHFP | `rdFingerprintGenerator` |
| 相似性 | Tanimoto、Dice、Cosine、Tversky 等 | `DataStructs` |
| 反应 | 反应 SMARTS、反应执行、产物枚举、反应指纹 | `rdChemReactions` |
| 2D/3D | 坐标、构象、力场、对齐、形状和绘图 | `AllChem`、`Draw` |
| 药效团 | 特征定义、2D/3D 药效团、特征映射 | `ChemicalFeatures`、`Pharm2D` |
| 片段分析 | BRICS、RECAP、Murcko、Fraggle、MMPA、MCS | `BRICS`、`Scaffolds` |
| 数据分析 | Pandas、数据库、向量、聚类、评分和 ML 辅助 | `PandasTools`、`Dbase`、`ML` |

### 适用范围

RDKit 适合单分子分析、化合物库处理、药物化学筛选、结构去重、SAR 分析、反应模板执行、分子机器学习前处理和基础三维结构分析。它不等同于量子化学程序、分子动力学程序或实验谱图解析软件。

安装、版本和顶层模块见 [手册第 0–2 节](doc/RDKit_manual.md#0-如何阅读本文档)。

---

## 二、安装与导入

### 安装要点

- Conda 环境优先使用 `conda-forge`；pip 安装时确认解释器和运行环境一致。
- 安装后检查 `rdkit.__version__`，不同版本可能改变默认参数、模块位置和弃用别名。
- InChI、Avalon、Cairo、FreeSASA 等功能可能取决于构建选项或可选组件。
- RDKit 自带的 FDef、盐规则、模板等数据应通过 `RDConfig` 或 `RDPaths` 定位，不要硬编码安装路径。

### 推荐入口

| 任务 | 入口 |
|---|---|
| 日常分子操作 | `from rdkit import Chem` |
| 构象、力场、指纹和综合功能 | `from rdkit.Chem import AllChem` |
| 描述符 | `Descriptors`、`rdMolDescriptors` |
| 指纹与相似度 | `rdFingerprintGenerator`、`DataStructs` |
| 绘图 | `Chem.Draw` |
| 标准化 | `MolStandardize.rdMolStandardize` |
| 参数和底层能力 | 对应 `rd*` 模块 |

`AllChem` 是历史兼容和便捷聚合入口；按具体模块导入有利于明确依赖和查找 API。

---

## 三、核心分子对象

### `Mol`

`Mol` 表示一个完整分子图，保存原子、键、环、立体、属性和构象。常用能力是获取原子/键/重原子/构象数量，按索引访问原子和键，遍历邻接关系，获取环信息，匹配子结构，管理构象，读写分子属性，以及合并、拆分或提取结构。

### `Atom`

`Atom` 保存元素、原子序数、同位素、形式电荷、芳香性、杂化、价态、隐式氢、自由基、手性标签、PDB 残基信息和自定义属性。可以读取或修改元素、同位素、电荷、芳香性、手性、显式/隐式价态和原子属性。

### `Bond`

`Bond` 保存起止原子、键级、键型、方向、顺反立体、芳香性、共轭性、环成员和自定义属性。支持单、双、三、芳香、配位、离子和其它特殊键型。

### `Conformer`

`Conformer` 保存一套二维或三维坐标；同一分子可拥有多个构象。它支持单点坐标读写、坐标矩阵、构象编号、二维/三维标记以及构象增删。

### 环、SGroup 和属性

`RingInfo` 可得到 SSSR、原子环、键环、环成员数和指定大小的环。`SubstanceGroup` 表示聚合物、数据组和特殊结构附加信息。`PropertyMol` 用于保存分子级属性。

核心对象、方法和异常见 [手册第 3 节](doc/RDKit_manual.md#3-核心数据模型rdkitchemrdchem)。

---

## 四、结构输入、输出与格式转换

### 字符串格式

| 格式 | 作用 | 代表 API |
|---|---|---|
| SMILES | 确定分子的紧凑表示 | `MolFromSmiles`、`MolToSmiles` |
| SMARTS | 带条件的查询模板 | `MolFromSmarts`、`MolToSmarts` |
| InChI | 标准结构标识 | `MolFromInchi`、`MolToInchi` |
| InChIKey | 外部数据库索引 | `MolToInchiKey` |
| CXSMILES | 坐标、映射和扩展属性 | `MolToCXSmiles` |

### 文件格式和供应器

`rdmolfiles` 支持 MOL/MolBlock、SDF、PDB、CML、CDXML、PNG 元数据和多种字符串供应器。`SDMolSupplier`、`ForwardSDMolSupplier`、`SmilesMolSupplier` 和 PDB 供应器适合批量读取；`SDWriter`、`SmilesWriter` 和 `PDBWriter` 负责写出。

### 输入输出注意事项

- 解析后通常需要 sanitize；异常结构可先关闭 sanitize，再手动处理。
- 批量供应器可能返回 `None`，必须保留坏记录和原始文本。
- 写出时应明确是否保留立体、同位素、显式氢、原子映射、坐标和属性。
- PDB 的键级可能不完整，可用 `rdDetermineBonds` 根据坐标推断连接和键级。
- JSON、SDF 和 PNG 元数据适合跨进程或带属性传输，但应检查版本兼容。

详见 [手册第 4 节](doc/RDKit_manual.md#4-io-模块rdmolopsrdmolfilesrddepictorrdcoordgen) 和 [第 5 节](doc/RDKit_manual.md#5-smiles--smarts--inchi--cxsmiles-全套)。

---

## 五、分子清洗、编辑与图操作

### 清洗和派生状态

`rdmolops` 负责使分子进入可计算状态：

- `SanitizeMol`：检查价态、芳香性、共轭、环和隐式氢。
- `Kekulize`、`SetAromaticity`：转换或重建芳香表示。
- `SetConjugation`、`SetHybridization`、`SetImplicitHs`：刷新派生属性。
- `AssignStereochemistry`：分配立体和 CIP 信息。
- `Cleanup`、`CleanupChirality`、`CleanupOrganometallics`：处理常见异常。

### 氢、片段和编辑

`AddHs`、`RemoveHs` 管理显式/隐式氢；`RWMol`、`EditableMol` 支持增删原子和键；`CombineMols` 合并分子；`GetMolFrags`、`SplitMolByPDBResidues` 拆分片段；`DeleteSubstructs`、`ReplaceSubstructs`、`FragmentOnBonds` 进行基于查询或键的编辑。

### 环、立体和查询调整

RDKit 可计算 SSSR 和对称 SSSR，寻找潜在立体中心，恢复楔形键、双键顺反和 3D 手性。`AdjustQueryProperties` 及相关函数可调整查询中的环、链、虚拟原子、连接点、度数和原子映射约束。

详细操作见 [手册第 3、4、6 节](doc/RDKit_manual.md#3-核心数据模型rdkitchemrdchem)。

---

## 六、SMILES、SMARTS 与子结构搜索

### SMILES

RDKit 支持规范/非规范 SMILES、带立体 SMILES、Kekulé SMILES、无立体 SMILES、同位素和 CXSMILES。输出时可控制规范化、立体、同位素、显式氢、原子映射、坐标和扩展字段。

### SMARTS

SMARTS 支持元素、原子序数、同位素、芳香性、杂化、价态、电荷、氢数、环、环大小、连接度、键级、方向、共轭、原子映射、递归查询和逻辑组合。

### 匹配 API

- `HasSubstructMatch`：判断是否存在匹配。
- `GetSubstructMatch`：获取首个原子映射。
- `GetSubstructMatches`：获取全部映射。
- `rdqueries`：构造类型化查询原子、查询键、比较条件和 `And/Or/Xor/Not` 逻辑。
- `SubstructLibrary`：在大型分子集合中进行指纹预筛选和子结构检索。
- `SynthonSpaceSearch`：面向反应构建块和合成空间搜索。

完整查询规则见 [手册第 5、6、14 节](doc/RDKit_manual.md#6-子结构查询-rdqueries-与-smarts-用法)。

---

## 七、化学反应与结构变换

`rdChemReactions` 用反应 SMARTS 表示反应物、试剂和产物的原子映射关系。主要功能包括创建、读取、验证、标准化和输出反应模板，执行反应物组合，枚举多组匹配产物，处理原子映射，以及计算反应差异指纹和结构指纹。

典型入口是 `ReactionFromSmarts`、反应对象的模板访问方法和 `RunReactants`。产物可能有多个组合，也可能需要重新 sanitize、去重和检查价态。

反应 SMARTS 描述的是图变换，不自动保证实验条件、产率、区域选择性、立体选择性和真实可行性。详细反应能力见 [手册第 7 节](doc/RDKit_manual.md#7-化学反应-rdchemreactionsallchemreaction) 和 [第 34 节](doc/RDKit_manual.md#34-高级示例-流水线-pipeline-与定制反应)。

---

## 八、分子描述符与性质计算

### 描述符类别

| 类别 | 主要指标 |
|---|---|
| 质量 | 分子量、精确质量、单同位素质量、分子式 |
| 极性 | TPSA、Labute ASA、VSA 分区 |
| 脂溶性 | Crippen LogP、MolMR、SlogP-VSA |
| 氢键 | HBD、HBA、Lipinski HBD/HBA |
| 组成 | 重原子、杂原子、价电子、形式电荷 |
| 柔性 | 可旋转键、环、桥环、螺环 |
| 立体 | 手性中心、已定义/未定义立体、CSP3 |
| 芳香性 | 芳香环、芳香碳环、芳香杂环 |
| 拓扑 | BertzCT、BalabanJ、Chi、Kappa、Hall-Kier |
| 片段 | 醇、胺、羰基、卤素、酰胺、醚和其它基团计数 |

### 主要模块

- `Descriptors`：常用描述符的统一 Python 入口。
- `rdMolDescriptors`：高性能 Boost 描述符、3D 指标和部分指纹。
- `MolSurf`：TPSA、LabuteASA 和表面积分区。
- `Crippen`：LogP 与分子摩尔折射率。
- `Lipinski`：Lipinski 规则相关计数。
- `QED`：定性估计药物相似性评分及其分量。
- `GraphDescriptors`：拓扑复杂度和图指数。
- `Fragments`：基于预定义 SMARTS 的功能团计数。

描述符是规则和模型指标，不等同于实验值；盐处理、氢键定义、芳香性模型和版本差异会影响结果。完整清单见 [手册第 8 节](doc/RDKit_manual.md#8-描述符-descriptorsrdmoldescriptorsmolsurfcrippenlipinskiqedgraphdescriptorsfragments)。

---

## 九、指纹、相似度和多样性

### 指纹类型

| 指纹 | 特点和用途 |
|---|---|
| Morgan/ECFP | 原子环境半径展开，适合相似性和 QSAR |
| RDKit Fingerprint | 路径和子结构编码，适合通用检索 |
| MACCS | 固定结构键集合，易解释和兼容 |
| AtomPair | 原子类型与拓扑距离 |
| Topological Torsion | 原子序列和拓扑扭转 |
| Avalon | Avalon 结构模式编码 |
| MHFP | MinHash 片段集合，适合近似检索 |
| Pharm2D | 药效团特征与距离 |
| EState | 电拓扑状态 |
| Reduced Graph/ErG | 缩减图模式 |
| 反应指纹 | 反应前后结构变化 |

### 生成和解释

`rdFingerprintGenerator` 以生成器对象统一配置半径、位数、计数模式、原子和键环境。Morgan、RDKit、AtomPair、TopologicalTorsion 和 Pattern 都有对应生成器。部分指纹可追踪 bit 对应的原子环境，适合解释命中原因；计数指纹保留频次，位向量适合快速集合运算。

### 相似度

`DataStructs` 支持显式位向量、稀疏向量和实值向量。主要指标是 Tanimoto、Dice、Cosine、Tversky、Sokal、Russel、Kulczynski、McConnaughey、Rogot-Goldberg 和 Braun-Blanquet；`Bulk*Similarity` 适合批量比较，邻居接口适合阈值检索。

### 多样性筛选

`SimDivFilters` 提供 `MaxMinPicker`、`LeaderPicker`、`HierarchicalClusterPicker` 和多种 `ClusterMethod`。它们可从大分子库中选择代表集、按相似度阈值分群或在聚类后选样。

详见 [手册第 9、26、28 节](doc/RDKit_manual.md#9-分子指纹与相似性-allchem--rdfingerprintgenerator--maccskeys--datastructs)。

---

## 十、2D 坐标和分子绘图

`rdDepictor` 和 `rdCoordGen` 根据拓扑生成二维坐标，支持规范朝向、模板匹配、坐标清理、结构拉直和系列分子共同骨架布局。

`Chem.Draw` 支持单分子图、网格图、反应图、PNG、SVG、PIL 和 Jupyter/IPython 展示，也支持原子、键、子结构、索引、立体标记、注释和属性高亮。

`rdMolDraw2D` 提供 `MolDraw2D`、`MolDraw2DSVG`、`MolDraw2DPNG`、`MolDraw2DCairo` 和 `MolDrawOptions`，可控制线宽、字体、颜色、画布、原子调色板、注释、背景和偏移。PNG 元数据接口可在图像中嵌入结构信息。

绘图和 2D 坐标见 [手册第 4、10、27 节](doc/RDKit_manual.md#27-分子绘图-rdkitchemdraw-与-rdchemdraw)。

---

## 十一、3D 构象、几何和力场

### 构象嵌入

`rdDistGeom` 和 `AllChem` 使用距离几何、ETKDG、ETKDGv2、ETKDGv3 生成单个或多个构象。参数可以控制随机种子、环约束、手性约束、距离约束、失败重试、坐标初始化和嵌入策略。必须检查嵌入返回状态和构象几何。

### 坐标与几何

`rdMolTransforms` 支持键长、键角、二面角、坐标点、平移、旋转和构象变换。`Geometry` 提供 `Point2D`、`Point3D`、`PointND`、规则网格、二面角、网格质心、突出距离、Tanimoto 和 Tversky 几何比较。`DistanceGeometry` 处理边界矩阵、三角平滑和距离约束。

### UFF 和 MMFF

`rdForceFieldHelpers` 和 `AllChem` 支持 UFF、MMFF94、MMFF94s 的参数检查、力场创建、能量计算、单构象优化、多构象优化和固定/位置/距离/角度/扭转约束。通用 `ForceField` 提供能量、迭代优化、约束和点数接口。

力场结果依赖参数、构象和力场模型；金属、自由基、特殊价态或参数缺失可能失败。力场优化不等同于量子化学优化。

详细说明见 [手册第 10、11、23、24 节](doc/RDKit_manual.md#10-三维构象-嵌入优化坐标生成)。

---

## 十二、3D 对齐、形状和药效团

### 分子对齐

`rdMolAlign` 支持基于原子映射、公共子结构或模板的三维对齐、RMSD、最佳构象搜索、变换矩阵和批量对齐。`rdMolTransforms` 负责对齐后的坐标变换和局部几何操作。

### 形状分析

`rdShapeAlign`、`rdShapeHelpers`、`rdGaussianShape` 和 `Subshape` 可计算三维体积、形状距离、高斯密度重叠、形状查询和形状叠合。所有结果都依赖有效构象、坐标尺度和构象选择。

### 化学特征和药效团

`ChemicalFeatures` 通过 FDef 识别氢键供体/受体、疏水中心、芳香环、正负电荷、可旋转键和自定义 SMARTS 特征。`Pharm2D` 将特征及拓扑距离编码为二维药效团指纹；`Pharm3D` 使用三维距离和方向；`FeatMaps` 进行特征匹配、映射和评分。

详见 [手册第 12、13 节](doc/RDKit_manual.md#12-药效团与化学特征-chemicalfeaturespharm2dpharm3dfeatmaps)。

---

## 十三、结构库、R 基团、骨架和片段

### 大型结构库

`SubstructLibrary` 结合分子、指纹和键索引，支持添加分子、快速匹配、命中计数、限制结果、序列化和检索顺序控制。`SynthonSpaceSearch` 面向构建块、合成空间和反应命中搜索。

### R 基团分解

`rdRGroupDecomposition` 根据一个或多个核心拆分取代基，支持原子/键比较策略、R 位点约束、核心匹配索引、行式或列式结果以及处理评分，适合 SAR 和系列化合物分析。

### 骨架和片段

- `MurckoScaffold`：提取环系统、连接体和骨架。
- `MakeScaffoldGeneric`：将骨架原子类型泛化。
- `rdScaffoldNetwork`：建立骨架网络。
- `BRICSDecompose`、`BRICSBuild`：可合成片段分解与重建。
- `RecapDecompose`：按 RECAP 规则建立层级片段树。
- `FragmentCatalog`：片段目录和片段匹配。
- `Fraggle`：可解释共享片段相似性。
- `rdMMPA`：匹配分子对和片段替换。

### 最大公共子结构

`rdFMCS.FindMCS` 支持原子、键、环和立体比较策略，可设置超时、最小结构、SMARTS 输出和匹配规则。`rdRascalMCES` 提供另一种 MCES 实现；历史 `Chem.MCS` 不建议新代码使用。

详见 [手册第 14、17、18、28 节](doc/RDKit_manual.md#14-子结构搜索-substructlibrarysynthonspacesearchrgroupdecomposition)。

---

## 十四、标准化、去盐、互变和过滤

### 标准化

`rdMolStandardize` 可组合完整预处理流程：`Cleanup` 清理异常，`Normalizer` 规范官能团，`Reionizer` 再离子化，`Uncharger` 中和电荷，`FragmentRemover` 删除指定片段，`LargestFragmentChooser` 保留最大片段，`DisconnectOrganometallics` 断开有机金属，`ValidateSmiles` 验证结构，`StandardizeSmiles` 一次性标准化，`Pipeline` 管理多阶段处理和状态日志。

### 盐、互变和查询

`SaltRemover` 按规则剥离常见盐和小片段。`TautomerEnumerator` 枚举并评分互变异构体，`TautomerQuery` 将互变关系作为查询模板。互变处理可能改变氢位置、电荷、键级和规范 SMILES，数据库去重前应固定策略。

### 结构过滤

`FilterCatalog` 支持 PAINS、BRENK、NIH、ZINC、ChEMBL、SureChEMBL 等结构警示目录，能够报告命中规则、描述和匹配细节。过滤结果是结构风险提示，不是毒性或活性结论。

标准化和过滤见 [手册第 15、16 节](doc/RDKit_manual.md#15-过滤与分子质量-filtercatalog-painsbrenknihlipinskiqed)。

---

## 十五、立体、哈希和结构身份

### 立体与共振

`EnumerateStereoisomers` 和 `rdEnumerateStereoisomers` 可枚举未定义手性中心、双键立体和环立体。`StereoEnumerationOptions` 控制唯一性、未定义立体、嵌入和枚举限制。RDKit 还支持 CIP、增强立体组、多种芳香性模型和共振标志。

### 哈希、键和序列化

`RegistrationHash` 可根据规范 SMILES、互变、立体、分子式、SGroup 等层生成注册哈希；`MolKey` 提供另一组分子身份层；`rdMolHash` 提供原子和结构层哈希；`rdMolInterchange` 支持 JSON 与分子对象互换。

这些工具适合去重、注册号、分层身份和数据交换，但应同时保存原始结构、标准化结构、版本和哈希策略。

详见 [手册第 19、21 节](doc/RDKit_manual.md#19-谱图与哈希-registrationhashmolkey原子指纹)。

---

## 十六、3D 描述符、网格和物理化学

### 3D 描述符

`Descriptors3D` 和 `rdMolDescriptors` 提供 PMI1/2/3、NPR1/2、Asphericity、Eccentricity、SpherocityIndex、RadiusOfGyration、InertialShapeFactor、PBF、USR、USRCAT、MORSE、WHIM、GETAWAY、RDF、Coulomb Matrix 和 EEM 等指标。它们都需要有效三维构象，部分指标对构象选择非常敏感。

### 网格、FreeSASA 和 EHT

`UniformGrid3D` 与 `UniformRealValueGrid3D` 表示三维网格，支持质心、端点、突出距离、Tanimoto/Tversky 比较和文件输出。`rdFreeSASA.CalcSASA` 计算溶剂可及表面积，支持原子半径、SASA 算法和极性分类。`rdEHTTools` 提供扩展 Hückel 结果，如轨道、原子电荷、键级、分子轨道和能量。

3D 描述符和物理化学模块是定性或快速计算工具，不替代高等级量子化学。详见 [手册第 20、23、29 节](doc/RDKit_manual.md#20-字符距离--训练-测试-切分--化学信息距离-descriptors3dspacialscore)。

---

## 十七、图论、向量和机器学习

### 图论和缩减图

`Graphs` 提供邻接矩阵和原子距离；`rdMolProcessing` 支持分子合并、子结构删除、子结构提取和氢处理；`ReducedGraphs` 将分子压缩为伪原子/伪键组成的图，并可生成 ErG 指纹。

### `DataStructs`

提供显式位向量、稀疏位向量、整数/浮点稀疏向量、离散向量、位运算、相似度、批量相似度、扁平距离矩阵、NumPy 转换和相似度邻居搜索。

### `rdkit.ML`

子系统包含 `Cluster` 聚类、`Data` 数据矩阵、`Descriptors` 统计描述符、`InfoTheory` 熵和信息量、`MatOps` 矩阵操作、`Scoring` AUC 等评分、`SLT` 决策树/统计学习和 `MLUtils` 数据切分辅助。它是基础工具层，不替代现代机器学习框架。

详见 [手册第 22、25、26 节](doc/RDKit_manual.md#25-ml-与统计学工具-rdkitml-子系统)。

---

## 十八、Pandas、数据库和批量处理

### Pandas

`PandasTools` 支持从 SMILES 列生成 `ROMol` 列，在 DataFrame 中渲染结构，导出 SDF，生成图像字符串，处理反应列，以及读取 SDF 为表格。`PandasPatcher` 提供 pandas 行为扩展。

### 数据库

`Dbase` 提供数据库连接、游标、数据库信息和 SQL 批量操作；`MolDb`、`MolCatalog` 管理分子条目和相似性数据。数据库建议同时保存原始输入、标准化结构、规范键、哈希、版本和错误记录。

### 批量建议

先解析并过滤无效结构，再统一盐、互变、电荷和氢策略；随后分批计算描述符和指纹，保留输入行号、异常原因和运行参数，最后输出表格、SDF、PNG 或数据库记录。

详见 [手册第 30、31 节](doc/RDKit_manual.md#30-pandas-集成-pandastools--pandaspatcher)。

---

## 十九、Avalon、MHFP、MIF、MMPA 和坐标推键

| 模块 | 功能 |
|---|---|
| `Avalon.pyAvalonTools` | Avalon 位指纹、计数指纹、位操作和相关指纹 |
| `rdMHFPFingerprint` | MinHash 分子片段指纹 |
| `rdMIF` | Molecular Interaction Field 分子相互作用场 |
| `rdMMPA` | Matched Molecular Pair 分子对和片段变化 |
| `rdDetermineBonds` | 从坐标推断连接关系和键级 |
| `rdtrajectory` | 分子轨迹和多帧结构 |
| `TemplateAlign` | 基于模板的对齐 |
| `SATIS` | 结构原子类型和拓扑分类 |
| `FragmentMatcher` | 片段匹配 |
| `EnumerateHeterocycles` | 杂环枚举 |

这些功能属于专项工具，参数和模块可用性应以当前安装版本为准。详见 [手册第 28、36 节](doc/RDKit_manual.md#28-mlavalonmhfpmif)。

---

## 二十、工具、配置和专项模块

### 配置、路径、日志和随机数

`RDConfig`、`RDPaths` 提供 RDKit 数据、源码、库和 Contrib 路径；`RDLogger` 控制日志等级和警告；`RDRandom` 提供相关随机数包装。构象、枚举和抽样任务应记录随机种子。

### VLib 和辅助工具

`VLib` 是节点、供应器、过滤器和变换器组成的数据流框架；`utils`、`ChemUtils`、`TestRunner`、`conftest` 主要服务旧式工具、测试和示例。日常业务优先使用 `Chem`、`AllChem` 和具体 `rd*` 模块。

### 其它专项模块

完整安装包还包含 `Suppliers`、`PropertyMol`、`SimpleEnum`、`ShowMols`、`PyMol`、`IFG`、`SATIS`、`FragmentMatcher`、`MolCatalog`、`MolKey`、`Features`、`FeatMaps`、`Pharm2D`、`Pharm3D`、`EState`、`AtomPairs`、`Recap` 和多个兼容/历史模块。

模块地图和工具说明见 [手册第 2、32、33、36 节](doc/RDKit_manual.md#36-附录-bapi-索引按模块完整列出)。

---

## 二十一、主要模块速查

| 模块 | 主要功能 | 代表入口 |
|---|---|---|
| `Chem` | 顶层聚合入口 | `MolFromSmiles` |
| `rdchem` | 核心对象 | `Mol`、`Atom`、`Bond` |
| `rdmolops` | 清洗、编辑和图操作 | `SanitizeMol`、`AddHs` |
| `rdmolfiles` | 结构格式读写 | `MolFromMolBlock`、供应器 |
| `AllChem` | 综合便捷入口 | `EmbedMolecule` |
| `Descriptors` | 常用性质 | `MolWt`、`MolLogP` |
| `rdMolDescriptors` | 高性能描述符 | `CalcTPSA`、`CalcNumHBA` |
| `rdFingerprintGenerator` | 可配置指纹 | `GetMorganGenerator` |
| `DataStructs` | 向量和相似度 | `TanimotoSimilarity` |
| `rdChemReactions` | 反应模板和执行 | `ReactionFromSmarts` |
| `rdFMCS` | 最大公共子结构 | `FindMCS` |
| `rdRGroupDecomposition` | R 基团分解 | `RGroupDecomposition` |
| `rdSubstructLibrary` | 大规模子结构库 | `SubstructLibrary` |
| `MolStandardize` | 标准化和互变 | `Pipeline` |
| `FilterCatalog` | 结构警示过滤 | `FilterCatalog` |
| `Scaffolds` | Murcko 和骨架网络 | `GetScaffoldForMol` |
| `BRICS` / `Recap` | 片段分解和重建 | `BRICSDecompose` |
| `Draw` | 2D 图像和反应绘图 | `MolToImage` |
| `rdDistGeom` | 3D 嵌入 | `ETKDGv3` |
| `rdForceFieldHelpers` | UFF/MMFF | `MMFFOptimizeMolecule` |
| `rdMolAlign` | 3D 对齐 | `GetBestRMS` |
| `ChemicalFeatures` | 化学特征 | `BuildFeatureFactory` |
| `Descriptors3D` | 3D 形状指标 | `Asphericity` |
| `rdFreeSASA` | 溶剂可及表面积 | `CalcSASA` |
| `rdEHTTools` | 扩展 Hückel | `RunMol` |
| `Geometry` | 点和网格 | `Point3D` |
| `SimDivFilters` | 多样性选样 | `MaxMinPicker` |
| `ML` | 聚类、熵、矩阵和评分 | `CalcAUC` |
| `PandasTools` | DataFrame 集成 | `AddMoleculeColumnToFrame` |
| `Dbase` / `MolDb` | 数据库和分子目录 | `DbConnect` |

完整模块索引见 [手册附录 B](doc/RDKit_manual.md#36-附录-bapi-索引按模块完整列出)。

---

## 二十二、按任务选择入口

| 任务 | 推荐入口 |
|---|---|
| SMILES 解析和规范化 | `Chem.MolFromSmiles`、`Chem.MolToSmiles` |
| SDF 批量读取 | `SDMolSupplier` |
| 结构合法性 | `SanitizeMol`、`ValidateSmiles` |
| 去盐和主片段选择 | `SaltRemover`、`FragmentRemover`、`LargestFragmentChooser` |
| 电荷标准化 | `Uncharger`、`Reionizer` |
| 互变体处理 | `TautomerEnumerator`、`TautomerQuery` |
| SMARTS 搜索 | `MolFromSmarts`、`HasSubstructMatch` |
| 大型子结构库 | `SubstructLibrary` |
| 反应执行 | `ReactionFromSmarts`、`RunReactants` |
| 分子量/LogP/TPSA | `Descriptors`、`rdMolDescriptors` |
| HBD/HBA/Lipinski | `Lipinski`、`CalcNumHBD`、`CalcNumHBA` |
| QED | `QED.qed` |
| Morgan 指纹 | `GetMorganGenerator` |
| 相似度 | `TanimotoSimilarity`、`BulkTanimotoSimilarity` |
| 多样性筛选 | `MaxMinPicker`、`LeaderPicker` |
| 2D 绘图 | `Draw.MolToImage`、`MolsToGridImage` |
| 3D 构象 | `rdDistGeom`、`EmbedMolecule` |
| 力场优化 | `MMFFOptimizeMolecule`、`UFFOptimizeMolecule` |
| 3D 对齐 | `rdMolAlign` |
| 骨架和片段 | `MurckoScaffold`、`BRICS`、`Recap` |
| 公共结构 | `rdFMCS.FindMCS` |
| 立体枚举 | `EnumerateStereoisomers` |
| 3D 形状 | `Descriptors3D`、`rdGaussianShape` |
| SASA | `rdFreeSASA.CalcSASA` |
| DataFrame 分析 | `PandasTools` |

精确签名和参数见 [完整 API 索引](doc/RDKit_manual.md#36-附录-bapi-索引按模块完整列出)。

---

## 二十三、推荐工作流

### 结构标准化

读取结构 → 检查解析结果 → sanitize → 去盐或选择主片段 → 规范电荷和官能团 → 规范互变体 → 生成规范 SMILES 或注册哈希 → 保存原始结构、处理结果和日志。

### 化合物库筛选

读取并过滤无效分子 → 统一标准化 → 计算质量、LogP、TPSA、HBD、HBA 和可旋转键 → 应用 PAINS/BRENK 或自定义 SMARTS → 生成指纹 → 相似度检索或多样性选样 → 输出结构、属性和过滤原因。

### SAR 分析

确定共同核心 → 使用 MCS 或 Murcko 提取骨架 → 进行 RGroup 分解 → 统计 R 位点和取代基 → 建立 DataFrame → 结合描述符、指纹、聚类和模型分析系列关系。

### 3D 分析

标准化 → 加氢 → ETKDG 生成多构象 → UFF/MMFF 优化 → 按能量或 RMSD 筛选 → 对齐共同骨架 → 计算 3D 描述符、药效团或形状分数。

### 反应处理

准备反应 SMARTS → 检查原子映射和价态 → 读取反应物 → 执行反应 → sanitize 产物 → 去重和过滤 → 输出产物结构、映射和反应指纹。

---

## 二十四、使用边界与常见问题

### 版本和兼容

RDKit 同时包含 C++ 包装、Python 便利函数和历史兼容别名。部分能力已迁移到 `rdMolStandardize`、`rdFingerprintGenerator` 和 `rdFMCS` 等新模块。应以当前安装版本的函数签名和完整手册为准。

### 结构与实验事实

RDKit 的价态、芳香性、立体、标准化、QED、PAINS 和 LogP 都是算法或规则结果。合法分子不一定是真实稳定化合物；反应结果不代表实验可行；过滤命中不等同于毒性结论。

### 2D、3D 和构象

SMILES 主要保存连接关系和立体标签，不是可靠的三维构象。3D 描述符、形状、SASA 和 EHT 必须使用有效三维坐标；二维坐标不能直接用于物理化学三维分析。

### 批量失败处理

应处理 `MolFrom*` 返回 `None`、sanitize/Kekulé/价态异常、力场参数缺失、嵌入失败、反应无匹配、供应器坏记录、缺少 conformer 和可选后端不可用等情况，不能简单忽略全部警告。

### 可复现性

建议记录 RDKit/Python 版本、输入格式、标准化参数、指纹参数、随机种子、力场类型、构象编号、过滤规则版本和异常记录。

详细验证报告见 [手册第 37 节](doc/RDKit_manual.md#37-验证报告与生成工具)。

---

## 二十五、枚举、异常和完整 API 索引

常用枚举包括 `BondType`、`BondDir`、`BondStereo`、`ChiralType`、`HybridizationType`、`AromaticityModel`、`ClusterMethod`、`ResonanceFlags`、`AtomCompare`、`BondCompare`、`RingCompare`、`HashLayer`、`HashScheme`、`PipelineStage` 和 `PipelineStatus`。

常见异常包括 SMILES 解析、价态、芳香性、Kekulé、sanitize、环和文件读写异常。批量任务应保存失败输入和异常原因。

完整枚举表见 [手册附录 A](doc/RDKit_manual.md#35-附录-a关键枚举与常量表)；按模块完整 API 见 [手册附录 B](doc/RDKit_manual.md#36-附录-bapi-索引按模块完整列出)。

---

## 二十六、参考入口

- [RDKit 官方网站](https://www.rdkit.org)
- [RDKit 官方文档](https://www.rdkit.org/docs/)
- [RDKit C++ API](https://www.rdkit.org/docs/cppapi/)
- [RDKit 完整函数手册](doc/RDKit_manual.md)
- [快速开始](doc/RDKit_manual.md#1-快速开始)
- [模块地图](doc/RDKit_manual.md#2-顶层包结构与子模块地图)
- [完整 API 索引](doc/RDKit_manual.md#36-附录-bapi-索引按模块完整列出)
- [验证报告](doc/RDKit_manual.md#37-验证报告与生成工具)
- [RDKit 官方仓库](https://github.com/rdkit/rdkit)

## 结论

RDKit 能覆盖从结构字符串解析、标准化、描述符和指纹计算，到分子库搜索、反应执行、2D/3D 建模、药效团分析、骨架/片段分析、Pandas 集成和基础机器学习前处理的完整本地工作流。

本文档用于快速回答“RDKit 能做什么、应该从哪里找 API、不同模块适合什么任务”。需要精确参数、完整方法、枚举值和可运行示例时，请继续查阅 [doc/RDKit_manual.md](doc/RDKit_manual.md)。

---

## 附：总结主题与完整手册对应关系

| 本总结主题 | 完整手册 |
|---|---|
| 概览、安装、模块地图 | [第 0–2 节](doc/RDKit_manual.md#0-如何阅读本文档) |
| 核心对象、清洗和图操作 | [第 3–4 节](doc/RDKit_manual.md#3-核心数据模型rdkitchemrdchem) |
| SMILES、SMARTS、InChI、子结构 | [第 5–6 节](doc/RDKit_manual.md#5-smiles--smarts--inchi--cxsmiles-全套) |
| 反应和描述符 | [第 7–8 节](doc/RDKit_manual.md#7-化学反应-rdchemreactionsallchemreaction) |
| 指纹、相似度和 3D | [第 9–11 节](doc/RDKit_manual.md#9-分子指纹与相似性-allchem--rdfingerprintgenerator--maccskeys--datastructs) |
| 药效团、对齐、结构库 | [第 12–14 节](doc/RDKit_manual.md#12-药效团与化学特征-chemicalfeaturespharm2dpharm3dfeatmaps) |
| 过滤、标准化、骨架、MCS | [第 15–18 节](doc/RDKit_manual.md#15-过滤与分子质量-filtercatalog-painsbrenknihlipinskiqed) |
| 哈希、3D 描述符、立体 | [第 19–21 节](doc/RDKit_manual.md#19-谱图与哈希-registrationhashmolkey原子指纹) |
| 图论、几何、ML、多样性 | [第 22–26 节](doc/RDKit_manual.md#22-图论算法-rdmolprocessinggraphsscaffolds) |
| 绘图、专项指纹、物理化学 | [第 27–29 节](doc/RDKit_manual.md#27-分子绘图-rdkitchemdraw-与-rdchemdraw) |
| Pandas、数据库、工具和完整索引 | [第 30–37 节](doc/RDKit_manual.md#30-pandas-集成-pandastools--pandaspatcher) |

> 部分阅读器对中文标题锚点的处理不同；若链接无法跳转，可在完整手册中按章节编号搜索。 
