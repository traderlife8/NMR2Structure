# PubChemPy 库功能总结

> 文档版本对照：
> - 官方文档：<https://docs.pubchempy.org/en/latest/>
> - 本机安装：`pubchempy 1.0.5`（位置：`C:\ProgramData\anaconda3\Lib\site-packages\pubchempy.py`）
> - 文档整理日期：2026-07-21

## 一、库概览

PubChemPy 是对 PubChem **PUG REST**（Power User Gateway REST）Web 服务的 Pythonic 封装，用于在 Python 中以面向对象的方式访问 PubChem 公开化学数据库。PubChem 包含上亿条化学结构及其生物、物理、毒理数据；PubChemPy 把这些远程接口抽象成易于使用的 Python 对象。

- **作用**：封装 PubChem PUG REST API，处理 URL 构造、HTTP 请求、JSON 解析、错误转换等细节。
- **支持 Python**：3.9 及以上。
- **源码仓库**：<https://github.com/mcs07/PubChemPy>
- **官方文档**：<https://docs.pubchempy.org/en/latest/>

### 核心特性

- 按 **name / SMILES / InChI / SDF / InChIKey / formula** 等标识符检索 Compound、Substance、Assay。
- 检索标准化后的 Compound 记录及其完整属性集。
- 在 **CID、SMILES、InChI、InChIKey、SDF、IUPAC Name** 等格式之间互转。
- 获取 PubChem 计算的 **分子描述符、XLogP、TPSA、复杂度、氢键供体/受体数** 等。
- 获取 **881 位的 PubChem CACTVS 结构指纹**（用于子结构/相似性搜索）。
- 生成 **2D 与 3D 坐标**，并获取 3D 体积、构象 RMSD、药效团特征等。
- 获取化合物所有 **同义词（IUPAC / 通用名 / 商品名 / 注册号）**。
- 一键下载 **XML / ASNT/B / JSON / SDF / PNG（图谱）** 文件。
- 与 **pandas DataFrame** 无缝集成，构造属性表。
- 完整覆盖 PubChem PUG REST，底层仍可调用 `request()` 发起自定义请求。

## 二、安装

```bash
# 推荐方式
pip install pubchempy

# 或者 conda
conda install -c conda-forge pubchempy

# 或者克隆源码安装
git clone https://github.com/mcs07/PubChemPy.git
cd PubChemPy
pip install -e .
```

依赖：`requests`（间接通过 `urllib`）、可选 `pandas`。

## 三、模块导出清单（`dir(pubchempy)` 去除下划线）

```
API_BASE, Assay, Atom, BadRequestError, Bond, BondType, Compound,
CompoundIdType, CoordinateType, ELEMENTS, MethodNotAllowedError, NotFoundError,
PROPERTY_MAP, ProjectCategory, PubChemHTTPError, PubChemPyDeprecationWarning,
PubChemPyError, ResponseParseError, ServerBusyError, ServerError, Substance,
TimeoutError, UnimplementedError, compounds_to_frame, download, get, get_aids,
get_all_sources, get_assays, get_cids, get_compounds, get_json, get_properties,
get_sdf, get_sids, get_substances, get_synonyms, request, substances_to_frame
```

> 说明：`HTTPError` 在本地命名空间中已可见，但官方 API 文档并未单独介绍；`annotations`/`certifi`/`functools`/`logging`/`os` 等是模块内部引用。

---

## 四、核心对象

### 1. `Compound`（标准化的化学结构记录）

PubChem Compound 数据库中每条结构对应唯一 CID，是经过标准化、去重后的化合物记录。

**构造方式**

| 方法 | 说明 |
|------|------|
| `Compound(record)` | 直接由 PUG REST 返回的 dict 构造（一般不直接调用） |
| `Compound.from_cid(cid, **kwargs)` | 给定 CID 抓取记录（首选方式） |

**常用属性**（部分需要再发起一次额外请求，已内部缓存）

| 属性 | 类型 | 说明 |
|------|------|------|
| `cid` | int | PubChem 化合物 ID |
| `elements` | list[str] | 原子对应的元素符号列表 |
| `atoms` | list[`Atom`] | 原子对象列表（含坐标） |
| `bonds` | list[`Bond`] | 键对象列表 |
| `charge` | int | 形式电荷 |
| `coordinate_type` | str | 2D / 3D 坐标类型 |
| `molecular_formula` | str | 分子式 |
| `molecular_weight` | float | 分子量（g/mol） |
| `smiles` | str | **绝对 SMILES**（含立体/同位素），1.0.5 起取代 `isomeric_smiles` |
| `connectivity_smiles` | str | **连通性 SMILES**（无立体/同位素），1.0.5 起取代 `canonical_smiles` |
| `inchi` | str | 标准 IUPAC InChI |
| `inchikey` | str | 27 位 InChIKey |
| `iupac_name` | str | 首选 IUPAC 系统名 |
| `xlogp` | float | 计算的 XLogP（辛醇-水分配系数） |
| `exact_mass` | float | 精确质量（Da） |
| `monoisotopic_mass` | float | 单同位素质量（Da） |
| `tpsa` | float | 拓扑极性表面积（Å²） |
| `complexity` | float | Bertz/Hendrickson/Ihlenfeldt 分子复杂度 |
| `h_bond_donor_count` | int | 氢键供体数 |
| `h_bond_acceptor_count` | int | 氢键受体数 |
| `rotatable_bond_count` | int | 可旋转单键数 |
| `heavy_atom_count` | int | 重原子数 |
| `isotope_atom_count` | int | 标记同位素原子数 |
| `atom_stereo_count` / `defined_atom_stereo_count` / `undefined_atom_stereo_count` | int | 原子立体中心计数（总计/已定义/未定义） |
| `bond_stereo_count` / `defined_bond_stereo_count` / `undefined_bond_stereo_count` | int | 键立体中心计数 |
| `covalent_unit_count` | int | 共价单元数 |
| `volume_3d` | float | 第一个多样构象的解析体积 |
| `conformer_rmsd_3d` | float | 构象采样 RMSD（Å） |
| `effective_rotor_count_3d` | int | 有效旋转键数 |
| `pharmacophore_features_3d` | list[str] | 3D 药效团特征列表 |
| `fingerprint` | str | 原始 hex 编码 881 位指纹 |
| `cactvs_fingerprint` | str | 881 位二进制结构指纹（按 CACTVS 子结构分段） |
| `synonyms` | list[str] | 排名后的所有同义词（需额外请求，结果缓存） |
| `sids` | list[int] | 关联的 Substance ID（额外请求，缓存） |
| `aids` | list[int] | 关联的 BioAssay ID（额外请求，缓存） |
| `record` | dict | PUG REST 返回的完整原始 dict |

**方法**

- `to_dict(properties=None)` → `dict`：导出属性 dict；`properties` 不指定时省略 `synonyms/aids/sids` 与已废弃的 `canonical_smiles/isomeric_smiles`。
- `to_series(properties=None)` → `pandas.Series`：同上，但返回 pandas Series。

### 2. `Atom`（原子对象）

```python
Atom(aid, number, x=None, y=None, z=None, charge=0)
```

- 属性：`aid`、`number`、`x`、`y`、`z`、`charge`
- `element`（property）→ 元素符号
- `coordinate_type`（property）→ `"2d"` 或 `"3d"`
- `set_coordinates(x, y, z=None)`
- `to_dict()`

### 3. `Bond`（键对象）

```python
Bond(aid1, aid2, order=BondType.SINGLE, style=None)
```

- 属性：`aid1`、`aid2`、`order`（`BondType` 枚举）、`style`
- `to_dict()`

### 4. `Substance`（原始沉积记录）

- 构造：`Substance.from_sid(sid, **kwargs)`
- 属性：
  - `sid`（int）、`synonyms`（list[str]）、`source_name`、`source_id`、`standardized_cid`
  - `standardized_compound`（property，对应 `Compound`，需额外请求、缓存）
  - `deposited_compound`（property，未经标准化的 `Compound`，无 CID）
  - `cids`（property，关联 CID 列表，额外请求）、`aids`（property，关联 AID 列表，额外请求）
- 方法：`to_dict(properties=None)`、`to_series(properties=None)`

### 5. `Assay`（生物测定记录）

- 构造：`Assay.from_aid(aid, **kwargs)`
- 属性：`aid`、`name`、`description`、`project_category`（`ProjectCategory` 枚举）、`comments`、`results`、`target`、`revision`、`aid_version`
- 方法：`to_dict(properties=None)`

---

## 五、检索函数（Search Functions）

### 1. `get_compounds(identifier, namespace='cid', searchtype=None, as_dataframe=False, **kwargs)`

- 作用：按给定标识符检索化合物记录。
- `identifier`：`str | int | list`
- `namespace`：`cid | name | smiles | sdf | inchi | inchikey | formula`
- `searchtype`：`substructure | superstructure | similarity | None`
- 返回：`list[Compound]` 或 `pd.DataFrame`（`as_dataframe=True`）
- 未命中返回空列表（自 v1.0.2 起）。

```python
import pubchempy as pcp
pcp.get_compounds("Aspirin", "name")                    # 按名字
pcp.get_compounds("CC(=O)Oc1ccccc1C(=O)O", "smiles")     # 按 SMILES
pcp.get_compounds("C9H8O4", "formula")                   # 按分子式
pcp.get_compounds("c1ccccc1", "smiles", searchtype="substructure")
pcp.get_compounds("aspirin", "name", as_dataframe=True)  # 直接得 DataFrame
```

### 2. `get_substances(identifier, namespace='sid', as_dataframe=False, **kwargs)`

- 检索 Substance 记录。
- `namespace`：`sid | name | sourceid/<source name>`
- 返回：`list[Substance]` 或 `pd.DataFrame`

### 3. `get_assays(identifier, namespace='aid', **kwargs)`

- 检索 BioAssay 记录。
- 返回：`list[Assay]`

### 4. `get_properties(properties, identifier, namespace='cid', searchtype=None, as_dataframe=False, **kwargs)`

- 按需拉取指定属性的轻量级检索，返回 `list[dict]` 或 `pd.DataFrame`，适合批量建表。
- `properties` 接受 `PROPERTY_MAP` 中所有键（含 `molecular_weight/xlogp/tpsa/smiles/inchi/iupac_name/...`）。

```python
pcp.get_properties(['MolecularFormula','MolecularWeight','XLogP','CanonicalSMILES'],
                   'aspirin', 'name')
```

### 5. `get_synonyms(identifier, namespace='cid', domain='compound', searchtype=None, **kwargs)`

- 拉取同义词列表（系统名/通用名/商品名/注册号）。
- 返回：`list[dict]`，每项含标识符与同义词数组。

---

## 六、标识符函数（Identifier Functions）

用于在 CID/SID/AID 与其它标识符之间互转。

### 1. `get_cids(identifier, namespace='name', domain='compound', searchtype=None, **kwargs)`
- 返回 `list[int]`，给定名字/SMILES/InChI 等得到对应 CID。
- `domain`：`compound | substance | assay`

### 2. `get_sids(identifier, namespace='cid', domain='compound', searchtype=None, **kwargs)`
- 返回 `list[int]`，由 CID/SMILES/name 等得到关联的 SID。

### 3. `get_aids(identifier, namespace='cid', domain='compound', searchtype=None, **kwargs)`
- 返回 `list[int]`，找到与该化合物/物质相关的所有 BioAssay ID。

### 4. `get_all_sources(domain='substance')`
- 返回 `list[str]`，列出当前所有 Substance 或 Assay 的提交方名称，可与 `get_substances(..., namespace='sourceid/<name>')` 配合使用。

---

## 七、请求函数（Request Functions）

底层 PUG REST 调用封装；高级用户可绕过 `Compound.from_cid` 等便捷方法，直接拿到原始响应。

### 1. `download(outformat, path, identifier, namespace='cid', domain='compound', operation=None, searchtype=None, overwrite=False, **kwargs)`
- 将结果直接保存为文件。
- `outformat`：`XML | ASNT/B | JSON | SDF | CSV | PNG | TXT`
- `path`：目标文件路径（`str` 或 `PathLike`）

### 2. `request(identifier, namespace='cid', domain='compound', operation=None, output='JSON', searchtype=None, **kwargs)`
- 返回 `http.client.HTTPResponse`，构造请求但不做高层解析。

### 3. `get(identifier, namespace='cid', domain='compound', operation=None, output='JSON', searchtype=None, **kwargs)`
- 返回 `bytes`，自动处理 PubChem 异步轮询（`Waiting` → 轮询到 `Ready`）。

### 4. `get_json(identifier, namespace='cid', domain='compound', operation=None, searchtype=None, **kwargs)`
- 返回 `dict | None`，自动解析 JSON；未命中静默返回 `None`。

### 5. `get_sdf(identifier, namespace='cid', domain='compound', operation=None, searchtype=None, **kwargs)`
- 返回 `str | None`，自动从响应中抽出 SDF 内容；未命中返回 `None`。

---

## 八、Pandas 集成函数

| 函数 | 说明 |
|------|------|
| `compounds_to_frame(compounds, properties=None)` | 由 `Compound` 列表/单个对象构造 DataFrame |
| `substances_to_frame(substances, properties=None)` | 由 `Substance` 列表/单个对象构造 DataFrame |
| 各搜索函数的 `as_dataframe=True` 参数 | 在检索阶段直接得到 DataFrame |

```python
from pubchempy import compounds_to_frame
cids = [2244, 1983, 3672]
df = compounds_to_frame(pcp.get_compounds(cids, 'cid'),
                         properties=['cid','iupac_name','molecular_weight','xlogp','tpsa'])
```

> 若已安装 pandas，可直接在 `get_compounds/get_substances/get_properties` 入口设置 `as_dataframe=True` 跳过分两步。

---

## 九、常量与枚举

### 1. `API_BASE`
```python
API_BASE = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'
```

### 2. `ELEMENTS`
- `dict[int, str]`，原子序数 → 元素符号。
- 包含 1~118 号标准元素，并额外提供 PubChem 用于非标准实体的特殊编码：
  - `252: 'Lp'`（lone pair）
  - `253: 'R'`（R 基）
  - `254/255: '*'`（dummy / unknown）

### 3. `PROPERTY_MAP`
- `dict[str, str]`，属性小写下划线名 → PubChem API 大写属性名。
- 共 43 个键，常见的有：
  - 结构：`smiles`, `connectivity_smiles`, `inchi`, `inchikey`, `iupac_name`
  - 分子量：`molecular_weight`, `exact_mass`, `monoisotopic_mass`
  - 描述符：`xlogp`, `tpsa`, `complexity`, `charge`, `heavy_atom_count`
  - 氢键与立体：`h_bond_donor_count`, `h_bond_acceptor_count`, `atom_stereo_count`, ...
  - 3D：`volume_3d`, `conformer_rmsd_3d`, `effective_rotor_count_3d`, `feature_*_count_3d`, ...

### 4. `CompoundIdType`（`IntEnum`）
| 名称 | 值 | 含义 |
|------|----|------|
| `DEPOSITED` | 0 | 原始沉积记录 |
| `STANDARDIZED` | 1 | 标准化后的形式 |
| `COMPONENT` | 2 | 标准化复合物的组件 |
| `NEUTRALIZED` | 3 | 中和形式 |
| `MIXTURE` | 4 | 混合物的组份 |
| `TAUTOMER` | 5 | 预测的互变异构体 |
| `IONIZED` | 6 | 预测的 pKa 电离形式 |
| `UNKNOWN` | 255 | 未知 |

### 5. `BondType`（`IntEnum`）
`SINGLE=1`, `DOUBLE=2`, `TRIPLE=3`, `QUADRUPLE=4`, `DATIVE=5`, `COMPLEX=6`, `IONIC=7`, `UNKNOWN=255`

### 6. `CoordinateType`（`IntEnum`）
- 类型维度：`TWO_D=1`, `THREE_D=2`, `SUBMITTED=3`, `EXPERIMENTAL=4`, `COMPUTED=5`, `STANDARDIZED=6`, `AUGMENTED=7`, `ALIGNED=8`, `COMPACT=9`
- 单位：`UNITS_ANGSTROMS=10`, `UNITS_NANOMETERS=11`, `UNITS_PIXEL=12`, `UNITS_POINTS=13`, `UNITS_STDBONDS=14`, `UNITS_UNKNOWN=255`

### 7. `ProjectCategory`（`IntEnum`）
- `MLSCN=1`, `MLPCN=2`, `MLSCN_AP=3`, `MLPCN_AP=4`, `JOURNAL_ARTICLE=5`, `ASSAY_VENDOR=6`, `LITERATURE_EXTRACTED=7`, `LITERATURE_AUTHOR=8`, `LITERATURE_PUBLISHER=9`, `RNAIGI=10`, `OTHER=255`

---

## 十、异常体系

```
PubChemPyError                # 所有异常的基类（继承自 Exception）
├── ResponseParseError        # 响应内容无法解析
└── PubChemHTTPError          # 通用 HTTP 错误（code, msg, details）
    ├── BadRequestError       # 400：请求格式错误
    ├── NotFoundError         # 404：找不到记录
    ├── MethodNotAllowedError # 405：不允许的方法/MIME
    ├── ServerError           # 500：服务器内部错误
    ├── UnimplementedError    # 501：未实现
    ├── ServerBusyError       # 503：服务器繁忙，应稍后重试
    └── TimeoutError          # 504：请求超时（可能因请求过宽）
```

警告类：`PubChemPyDeprecationWarning`（继承自 `Warning`），例如 `canonical_smiles` / `isomeric_smiles` 自 1.0.5 起被标记弃用。

> 行为变更：
> - 自 v1.0.2 起，搜索函数在无结果时返回 `[]` 而非抛出异常；
> - `Compound.from_cid(invalid)` 仍会抛出 `NotFoundError`；
> - 自 v1.0.3 起，`Compound.atoms` 与 `Compound.bonds` 返回 `Atom`/`Bond` 对象列表，而非 dict。

---

## 十一、典型使用示例

### 1. 抓取单个化合物

```python
import pubchempy as pcp

c = pcp.Compound.from_cid(2244)        # 阿司匹林
print(c.iupac_name, c.molecular_weight, c.xlogp, c.tpsa)
print(c.smiles, c.cactvs_fingerprint[:60], "…")
```

### 2. 按名称/SMILES/分子式搜索

```python
pcp.get_compounds("glucose", "name")
pcp.get_compounds("CC(=O)Oc1ccccc1C(=O)O", "smiles")
pcp.get_compounds("C9H8O4", "formula")
```

### 3. 子结构 / 相似性 / 超结构搜索

```python
# 含苯环的所有化合物（CID 列表）
cids = pcp.get_cids("c1ccccc1", "smiles", searchtype="substructure")
# 与咖啡因相似度 90% 以上
sim = pcp.get_compounds("CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                        "smiles", searchtype="similarity", Threshold=90)
```

### 4. 同义词 / 关联 SID、AID

```python
c = pcp.Compound.from_cid(2244)
print(c.synonyms[:5])   # 自动发一次额外请求
print(c.sids, c.aids)
```

### 5. 用 pandas 构造属性表

```python
cids = [2244, 1983, 3672, 3672]
df = pcp.get_properties(
    ["MolecularFormula", "MolecularWeight", "XLogP", "TPSA", "CanonicalSMILES"],
    cids, "cid",
)
print(df)
```

### 6. 下载文件（PNG/SDF/JSON/XML）

```python
# 2D 结构 PNG
pcp.download("PNG", "aspirin.png", 2244)
# 3D 构象 SDF
pcp.download("SDF", "aspirin_3d.sdf", 2244, operation="record", overwrite=True)
```

### 7. 自定义 PUG REST 请求

```python
data = pcp.get_json("aspirin", "name", domain="compound")
cids = pcp.get_cids("aspirin", "name")
```

### 8. 检索 Substance 与 Assay

```python
sources = pcp.get_all_sources(domain="substance")
subs = pcp.get_subsources("KEGG", "sourceid/KEGG")  # 旧 API 形式
subs = pcp.get_substances("lipid maps", "name")

assay = pcp.Assay.from_aid(1000)
print(assay.name, assay.project_category)
```

---

## 十二、高级使用

| 主题 | 说明 |
|------|------|
| **避免 TimeoutError** | 504 通常因请求范围过宽；减小返回规模（如使用 `as_dataframe=True` + `get_properties` 而非 `get_compounds`），或拆分多次请求。 |
| **Logging** | 模块使用标准 `logging`：`logging.getLogger('pubchempy')`。 |
| **代理访问** | 在调用前注入 `urllib.request` 的 `ProxyHandler`，或设置环境变量 `HTTP_PROXY`/`HTTPS_PROXY`。 |
| **自定义请求** | 直接使用 `request/get/get_json/get_sdf/download`；遵循 PUG REST URL 规范 `domain/namespace/identifier/operation/output`。 |
| **PUG REST 限制** | 单次最多返回约 100,000 条记录；超出需用 PUG VIEW 或批量提交。 |

---

## 十三、常用入口速查表

| 任务 | 推荐调用 |
|------|---------|
| 已知 CID 拿全部属性 | `Compound.from_cid(cid)` |
| 按名称/SMILES/InChI 找化合物 | `get_compounds(id, namespace)` |
| 只取若干字段 | `get_properties(props, id, namespace)` |
| 找 CID | `get_cids(id, namespace)` |
| 拿所有同义词 | `get_synonyms(id, namespace)` |
| 拿 3D 构象 SDF | `download('SDF', path, cid, operation='record')` 或 `get_sdf(cid, operation='record')` |
| 拿 2D 图 PNG | `download('PNG', path, cid)` |
| 转成 DataFrame | 任意检索函数加 `as_dataframe=True`，或 `compounds_to_frame(...)` |
| 自定义 PUG REST | `request()` / `get()` / `get_json()` / `get_sdf()` |

---

## 十四、引用与链接

- 官方文档：<https://docs.pubchempy.org/en/latest/>
- 用户指南：
  - [Introduction](https://docs.pubchempy.org/en/latest/guide/introduction.html)
  - [Installation](https://docs.pubchempy.org/en/latest/guide/install.html)
  - [Getting started](https://docs.pubchempy.org/en/latest/guide/gettingstarted.html)
  - [Searching](https://docs.pubchempy.org/en/latest/guide/searching.html)
  - [Compounds](https://docs.pubchempy.org/en/latest/guide/compound.html)
  - [Substances](https://docs.pubchempy.org/en/latest/guide/substance.html)
  - [Properties](https://docs.pubchempy.org/en/latest/guide/properties.html)
  - [pandas integration](https://docs.pubchempy.org/en/latest/guide/pandas.html)
  - [Download](https://docs.pubchempy.org/en/latest/guide/download.html)
  - [Advanced usage](https://docs.pubchempy.org/en/latest/guide/advanced.html)
  - [PUG REST](https://docs.pubchempy.org/en/latest/guide/pugrest.html)
- API 文档：<https://docs.pubchempy.org/en/latest/api.html>
- PubChem PUG REST 规范：<https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest>
- PubChemPy 源码：<https://github.com/mcs07/PubChemPy>