# 化学反应配平 API

用 **Python + FastAPI + SymPy** 实现的化学反应式自动配平服务，支持普通反应、
**离子反应**与**含电子的半反应**，并支持可选的 **系数比约束（ratios）**；
另提供 **酸性半反应自动补全接口**：只给一种反应物和一种生成物的骨架，
H2O、H^+、e^- 由接口作为带符号未知项自动求解并放到正确的一侧。

按元素守恒与电荷守恒构建**精确整数矩阵**（SymPy 有理数运算，无浮点误差），仅当：

1. 矩阵零空间维数**恰为 1**，且
2. 该唯一解可以整体取号为**全部严格正整数**

时才返回配平结果；系数约成最大公约数为 1 的最小正整数，并返回反应两侧各元素的
原子总数与带符号总电荷。不会为了配平而自动增删任何物种。

可选的 `ratios` 数组把「某两种物质的系数比为 a:b」作为**精确矩阵行**追加进
守恒矩阵联立求解：原本有多组独立解的反应可能被约束收敛到唯一全正解；
约束矛盾或约束后无全正解返回 `unbalanceable`，约束后仍多解返回
`multiple_solutions`；省略 `ratios` 时行为与之前完全一致。

## 支持与不支持

**支持**

- 元素符号（1–118 号真实元素，首字母大写、最多一个小写字母）
- 正整数下标（如 `H2O`、`Fe2(SO4)3`），不允许前导 0
- **最多四层**圆括号嵌套，如 `K4(Fe(CN)6)`、`Mg2(Al(Si(BO4)2)3)2`
- 重复元素累计计数（括号内外同名元素相加）
- **电荷后缀**：化学式末尾的 `^+`、`^-`、`^n+`、`^n-`
  （n 为无前导零的正整数，省略时按 1 计），如 `Na^+`、`SO4^2-`、`Fe^3+`；
  无后缀表示中性
- **电子**：仅写作 `e^-`（原子计数为空、电荷为 -1），可配平半反应

**不支持**（直接判为非法化学式并给出字符位置）

- 裸 `+`/`-` 电荷记号（`Na+`、`OH-`）——请改写为 `Na^+`、`OH^-`
- 电子的其他写法（`e-`、`e^+`、`e^2-`）——电子仅 `e^-`
- 点号 / 中点结晶水（`CuSO4.5H2O`、`CuSO4·5H2O`）——请把 `CuSO4` 与 `5H2O` 拆成独立物质
- 方括号、花括号（请改写为圆括号）

## 一键启动（单条命令）

需要本机已安装 Docker 与 Docker Compose v2：

```bash
docker compose up --build
```

该命令在构建镜像时自动完成依赖安装（`pip install -r requirements.txt`），
无需任何手动安装步骤；加 `-d` 可后台运行。停止服务：

```bash
docker compose down
```

### 必要配置

| 配置项 | 位置 | 默认值 | 说明 |
| --- | --- | --- | --- |
| `BALANCER_PORT` | `docker-compose.yml` 环境变量 / `.env` | `8000` | 宿主机映射端口，如端口冲突可 `BALANCER_PORT=9000 docker compose up --build` |
| 容器内监听端口 | `Dockerfile` CMD | `8000` | 与 `EXPOSE` 保持一致，一般无需修改 |

### API 地址与端口

- 服务地址：`http://localhost:8000`
- 配平接口：`POST http://localhost:8000/balance`
- 酸性半反应自动补全接口：`POST http://localhost:8000/balance/acidic-half`
- 健康检查：`GET http://localhost:8000/health`
- Swagger 交互文档：`http://localhost:8000/docs`
- OpenAPI Schema：`http://localhost:8000/openapi.json`

## 本地开发（不用 Docker）

```bash
pip install -r requirements.txt
python -m uvicorn app.main:app --host 0.0.0.0 --port 8000
# 运行测试
pip install pytest && python -m pytest -q
```

## API 说明

### `POST /balance`

请求体：

```json
{
  "reactants": ["Fe", "O2"],
  "products": ["Fe3O4"]
}
```

- 两个数组均至少含 1 个、至多 50 个化学式；
- **数组顺序即返回系数顺序**，不做重排、不去重；
- 可选 `ratios` 数组（至多 100 项）指定系数比约束，见下文「系数比约束」。

#### 系数比约束 `ratios`（可选）

每项约束用 `left`、`right` 指定**两种不同物质**（`side` 取 `reactants` 或
`products`，`index` 为该侧从 0 开始的下标），并用 `ratio: [a, b]`（两个
**正整数**）要求两者的配平系数比为 `a : b`：

```json
{
  "reactants": ["C", "O2", "CO"],
  "products": ["CO2"],
  "ratios": [
    {
      "left":  {"side": "reactants", "index": 0},
      "right": {"side": "reactants", "index": 2},
      "ratio": [1, 1]
    }
  ]
}
```

`C + O2 + CO -> CO2` 本身存在多组独立解；约束 `c(C) : c(CO) = 1 : 1` 后
收敛到唯一解（HTTP 200，`balanced: true`）：

```json
{
  "balanced": true,
  "coefficients": {"reactants": [2, 3, 2], "products": [4]},
  "equation": "2C + 3O2 + 2CO -> 4CO2",
  "atom_totals": {
    "reactants": {"C": 4, "O": 8},
    "products": {"C": 4, "O": 8}
  },
  "charge_totals": {"reactants": 0, "products": 0}
}
```

约束非法时返回 **HTTP 422**（FastAPI 标准校验错误结构，`detail` 中给出
具体原因）：

- `side` 不是 `reactants` / `products`；
- `index` 越界（小于 0 或不小于该侧物质数）；
- `left` 与 `right` 指向同一物质；
- `ratio` 不是恰好两个正整数（含 0、负数、小数、长度不为 2）。

约束本身合法但求解失败时仍走原有的 `reason` 语义（HTTP 200，
`balanced: false`）：约束矛盾或约束后不存在全正解 → `unbalanceable`；
约束后零空间维数仍大于 1 → `multiple_solutions`。

#### 成功（HTTP 200，`balanced: true`）

```json
{
  "balanced": true,
  "reactants": [
    {"formula": "Fe", "coefficient": 3},
    {"formula": "O2", "coefficient": 2}
  ],
  "products": [
    {"formula": "Fe3O4", "coefficient": 1}
  ],
  "coefficients": {"reactants": [3, 2], "products": [1]},
  "equation": "3Fe + 2O2 -> Fe3O4",
  "atom_totals": {
    "reactants": {"Fe": 3, "O": 4},
    "products": {"Fe": 3, "O": 4}
  },
  "charge_totals": {"reactants": 0, "products": 0}
}
```

`atom_totals` 给出配平后两侧**各元素的原子总数**（两侧必然相等）；
`charge_totals` 给出两侧**带符号的总电荷**（两侧必然相等，中性反应为 0）。

离子 / 半反应示例（`MnO4^- + H^+ + e^- -> Mn^2+ + H2O`）：

```json
{
  "balanced": true,
  "coefficients": {"reactants": [1, 8, 5], "products": [1, 4]},
  "equation": "MnO4^- + 8H^+ + 5e^- -> Mn^2+ + 4H2O",
  "atom_totals": {
    "reactants": {"H": 8, "Mn": 1, "O": 4},
    "products": {"H": 8, "Mn": 1, "O": 4}
  },
  "charge_totals": {"reactants": 2, "products": 2}
}
```

#### 无法配平（HTTP 200，`balanced: false`）

零空间为空（原子或电荷无法同时守恒），或唯一解含零/正负混合系数，
无法得到全部严格正系数：

```json
{
  "balanced": false,
  "reason": "unbalanceable",
  "message": "元素守恒方程无解，该反应无法配平"
}
```

#### 多组独立解（HTTP 200，`balanced: false`）

零空间维数大于 1，通常因为输入混入了两个以上独立反应：

```json
{
  "balanced": false,
  "reason": "multiple_solutions",
  "message": "零空间维数为 2（大于 1），存在多组独立的配平系数"
}
```

> 「无法配平」与「存在多组独立解」通过 `reason` 字段明确区分，
> 不会返回任意一组解来冒充唯一答案。

#### 非法化学式（HTTP 400，`balanced: false`）

一次性返回两侧**全部**非法化学式（含非法电荷写法），`position` 为化学式中
首个非法字符的 **0 起始下标**（未闭合的 `(` 指向该括号本身；电荷后缀不完整时
指向应有后续字符的位置）：

```json
{
  "balanced": false,
  "reason": "invalid_formula",
  "message": "存在非法化学式，错误位置见 errors",
  "errors": [
    {
      "side": "reactants",
      "index": 0,
      "formula": "SO4^02-",
      "position": 4,
      "message": "电荷数必须为正整数且不能以 0 开头"
    }
  ]
}
```

### `POST /balance/acidic-half`

酸性条件下的**半反应自动补全**：只需给出仅含**一种反应物和一种生成物**的
半反应骨架（化学式与电荷写法与 `/balance` 相同），接口把 **H2O、H^+、e^-**
作为**可带符号的未知项**与两种骨架物质（系数限定为正数）一起联立元素守恒与
电荷守恒方程精确求解。解出后：**负值**未知项移到反应物侧、**正值**放到
生成物侧、**零项删去**，全部系数约成互质整数。

请求体（两个数组都必须**恰好含 1 个**化学式，否则 HTTP 422）：

```json
{
  "reactants": ["MnO4^-"],
  "products": ["Mn^2+"]
}
```

骨架中**不得包含 H2O、H^+、e^-**（它们由接口自动补全）；按解析结果匹配，
`HOH`、`H1^+` 等等价写法同样会被拒绝（HTTP 400，`reason: "forbidden_species"`）。

#### 成功（HTTP 200，`balanced: true`）

返回结构与 `/balance` 完全一致：骨架物质排在各自侧最前，自动补全的
H2O、H^+、e^- 按此顺序追加到实际所属的一侧：

```json
{
  "balanced": true,
  "reactants": [
    {"formula": "MnO4^-", "coefficient": 1},
    {"formula": "H^+", "coefficient": 8},
    {"formula": "e^-", "coefficient": 5}
  ],
  "products": [
    {"formula": "Mn^2+", "coefficient": 1},
    {"formula": "H2O", "coefficient": 4}
  ],
  "coefficients": {"reactants": [1, 8, 5], "products": [1, 4]},
  "equation": "MnO4^- + 8H^+ + 5e^- -> Mn^2+ + 4H2O",
  "atom_totals": {
    "reactants": {"H": 8, "Mn": 1, "O": 4},
    "products": {"H": 8, "Mn": 1, "O": 4}
  },
  "charge_totals": {"reactants": 2, "products": 2}
}
```

系数为 0 的自动补全项会被删去，例如 `Fe^2+ -> Fe^3+` 的结果为
`Fe^2+ -> Fe^3+ + e^-`（不出现 H2O 与 H^+）；H2O 解为负值时移到反应物侧，
例如 `SO3^2- + H2O -> SO4^2- + 2H^+ + 2e^-`。

#### 失败（`balanced: false`）

- **HTTP 200，`reason: "unbalanceable"`**：守恒方程无解，或唯一解中骨架
  物质系数为 0 / 异号，不存在唯一全正骨架解（如 `KMnO4 -> K2MnO4`，
  K 元素无法守恒）；
- **HTTP 200，`reason: "multiple_solutions"`**：零空间维数大于 1，骨架
  不能唯一配平（如 `H2O2 -> O2`，骨架只含 H、O 时 H2O 与 H^+ 无法唯一确定）；
- **HTTP 200，`reason: "no_electron_transfer"`**：唯一解中电子系数为 0，
  该骨架不是氧化还原半反应（如 `H2CO3 -> HCO3^-` 只是酸碱解离）；
- **HTTP 400，`reason: "invalid_formula"`**：化学式非法，结构同 `/balance`；
- **HTTP 400，`reason: "forbidden_species"`**：骨架含 H2O、H^+、e^-：

```json
{
  "balanced": false,
  "reason": "forbidden_species",
  "message": "骨架中不得包含 H2O、H^+、e^-（由接口自动补全），详见 errors",
  "errors": [
    {
      "side": "reactants",
      "index": 0,
      "formula": "H2O",
      "message": "H2O 由接口在酸性条件下自动补全，不得出现在半反应骨架中"
    }
  ]
}
```

## 可调用示例

### curl

```bash
# 普通反应
curl -s -X POST http://localhost:8000/balance \
  -H 'Content-Type: application/json' \
  -d '{"reactants": ["K4(Fe(CN)6)", "KMnO4", "H2SO4"],
       "products": ["KHSO4", "Fe2(SO4)3", "MnSO4", "HNO3", "CO2", "H2O"]}'

# 酸性条件下的半反应（含电子）
curl -s -X POST http://localhost:8000/balance \
  -H 'Content-Type: application/json' \
  -d '{"reactants": ["MnO4^-", "H^+", "e^-"], "products": ["Mn^2+", "H2O"]}'

# 离子方程式
curl -s -X POST http://localhost:8000/balance \
  -H 'Content-Type: application/json' \
  -d '{"reactants": ["Zn", "Cu^2+"], "products": ["Zn^2+", "Cu"]}'

# 系数比约束：C + O2 + CO -> CO2 原本多解，约束 c(C):c(CO)=1:1 后收敛到唯一解
curl -s -X POST http://localhost:8000/balance \
  -H 'Content-Type: application/json' \
  -d '{"reactants": ["C", "O2", "CO"], "products": ["CO2"],
       "ratios": [{"left": {"side": "reactants", "index": 0},
                   "right": {"side": "reactants", "index": 2},
                   "ratio": [1, 1]}]}'

# 酸性半反应自动补全：只给骨架，H2O、H^+、e^- 由接口求解
curl -s -X POST http://localhost:8000/balance/acidic-half \
  -H 'Content-Type: application/json' \
  -d '{"reactants": ["MnO4^-"], "products": ["Mn^2+"]}'
```

### Python（仅标准库）

`examples/example_call.py` 演示普通配平、离子/半反应、系数比约束、无法配平、
多组解、非法化学式、非法比值约束（422）与酸性半反应自动补全等场景：

```bash
python examples/example_call.py
# 服务地址不同时：BASE_URL=http://localhost:9000 python examples/example_call.py
```

### Python（requests）

```python
import requests

resp = requests.post(
    "http://localhost:8000/balance",
    json={"reactants": ["H2O"], "products": ["H2", "O2"]},
)
print(resp.json())
# 2H2O -> 2H2 + O2

# 带系数比约束：要求 c(H2) : c(O2) = 2 : 1
resp = requests.post(
    "http://localhost:8000/balance",
    json={
        "reactants": ["H2O"],
        "products": ["H2", "O2"],
        "ratios": [
            {
                "left": {"side": "products", "index": 0},
                "right": {"side": "products", "index": 1},
                "ratio": [2, 1],
            }
        ],
    },
)
print(resp.json())
# 2H2O -> 2H2 + O2（约束与守恒一致，仍收敛到同一最小整数解）

# 酸性半反应自动补全：骨架 MnO4^- -> Mn^2+
resp = requests.post(
    "http://localhost:8000/balance/acidic-half",
    json={"reactants": ["MnO4^-"], "products": ["Mn^2+"]},
)
print(resp.json())
# MnO4^- + 8H^+ + 5e^- -> Mn^2+ + 4H2O
```

## 算法与判定规则

1. 用递归下降解析器校验并解析每个物种，得到 `{元素: 原子数}` 与带符号电荷
   （电子 `e^-` 的原子计数为空、电荷为 -1）；
2. 以物质为列构建矩阵 `A`：每种元素一行，反应物列取正、生成物列取负；
   额外加入一行**带符号电荷行**（同一符号约定），中性反应时该行为全零，
   不改变原有行为；
3. 若请求带 `ratios`，把每条系数比约束 `c_left : c_right = a : b` 转成
   精确整数行 `b·c_left − a·c_right = 0` 追加到矩阵末尾，与元素、电荷
   守恒行联立；
4. `sympy.Matrix(...).nullspace()` 在精确有理数域上求零空间：
   - 零空间维数 0 → `unbalanceable`（含比值约束互相矛盾的情形）；
   - 维数 > 1 → `multiple_solutions`（约束后仍多解时同样返回）；
5. 唯一基向量含零分量或正负混合 → `unbalanceable`（不存在全部严格正系数）；
6. 整体取号为正，按分母公倍数转整数，再除以全部系数的最大公约数，
   得到 gcd = 1 的最小正整数解；
7. 校验 `A·c = 0`（原子守恒、电荷守恒且满足全部比值约束）后返回系数、
   两侧原子总数与两侧带符号总电荷。

`/balance/acidic-half` 的求解流程类似：矩阵列固定为
`[骨架反应物, 骨架生成物, H2O, H^+, e^-]`（骨架反应物列取正，其余列取负，
即未知项暂定生成物侧），行取骨架元素与 H、O 的并集再加一行电荷；
零空间维数必须恰为 1，且两种骨架物质系数同号非零（整体取号为正）；
H2O、H^+、e^- 的带符号解为负则移到反应物侧、为正放到生成物侧、
为零则删去；电子系数为 0 判为 `no_electron_transfer`；
最终全部系数约成互质整数并同样校验 `A·c = 0`。

## 项目结构

```
.
├── app/
│   ├── main.py       # FastAPI 入口（/balance、/balance/acidic-half、/health）
│   ├── schemas.py    # Pydantic 请求/响应模型（含 ratios 与半反应骨架校验）
│   ├── balancer.py   # 元素+电荷守恒矩阵、比值约束行、零空间配平与酸性半反应补全
│   └── parser.py     # 化学式/电荷/电子递归下降解析（四层括号/位置报错）
├── examples/
│   └── example_call.py
├── tests/            # pytest 测试（81 个用例）
├── Dockerfile
├── docker-compose.yml
└── requirements.txt
```
