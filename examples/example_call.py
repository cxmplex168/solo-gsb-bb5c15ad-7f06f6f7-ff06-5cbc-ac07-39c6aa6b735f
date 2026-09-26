"""配平 API 调用示例（仅依赖 Python 标准库）。

用法::

    # 先启动服务
    docker compose up --build

    # 再运行示例（默认 http://localhost:8000，可用环境变量 BASE_URL 覆盖）
    python examples/example_call.py
"""

from __future__ import annotations

import json
import os
import urllib.error
import urllib.request

BASE_URL = os.environ.get("BASE_URL", "http://localhost:8000")


def call_balance(
    reactants: list[str], products: list[str], ratios: list[dict] | None = None
) -> dict:
    payload_dict: dict = {"reactants": reactants, "products": products}
    if ratios is not None:
        payload_dict["ratios"] = ratios
    return _post(f"{BASE_URL}/balance", payload_dict)


def call_acidic_half(reactant: str, product: str) -> dict:
    return _post(
        f"{BASE_URL}/balance/acidic-half",
        {"reactants": [reactant], "products": [product]},
    )


def _post(url: str, payload_dict: dict) -> dict:
    payload = json.dumps(payload_dict).encode("utf-8")
    request = urllib.request.Request(
        url,
        data=payload,
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(request) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        # 非法化学式 / 骨架含自动补全物种返回 400、非法比值约束返回 422，
        # 响应体仍是 JSON
        return json.loads(exc.read().decode("utf-8"))


def show(
    title: str,
    reactants: list[str],
    products: list[str],
    ratios: list[dict] | None = None,
) -> None:
    print(f"=== {title} ===")
    print(f"请求: reactants={reactants}, products={products}, ratios={ratios}")
    result = call_balance(reactants, products, ratios)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    print()


def show_acidic_half(title: str, reactant: str, product: str) -> None:
    print(f"=== {title} ===")
    print(f"请求: reactants=[{reactant!r}], products=[{product!r}]")
    result = call_acidic_half(reactant, product)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    print()


def main() -> None:
    # 1. 普通配平：Fe + O2 -> Fe3O4  =>  3Fe + 2O2 -> Fe3O4
    show("普通配平", ["Fe", "O2"], ["Fe3O4"])

    # 2. 酸性条件下的半反应：MnO4^- + 8H^+ + 5e^- -> Mn^2+ + 4H2O
    show("含电子的半反应", ["MnO4^-", "H^+", "e^-"], ["Mn^2+", "H2O"])

    # 3. 离子方程式：Zn + Cu^2+ -> Zn^2+ + Cu
    show("离子方程式", ["Zn", "Cu^2+"], ["Zn^2+", "Cu"])

    # 4. 四层括号 + 重复元素累计：K4(Fe(CN)6) 参与的反应
    show(
        "复杂括号",
        ["K4(Fe(CN)6)", "KMnO4", "H2SO4"],
        ["KHSO4", "Fe2(SO4)3", "MnSO4", "HNO3", "CO2", "H2O"],
    )

    # 5. 无法配平：H2 + O2 -> H2O + Fe
    show("无法配平", ["H2", "O2"], ["H2O", "Fe"])

    # 6. 电荷不守恒：Fe^2+ -> Fe^3+（缺少电子）
    show("电荷不守恒", ["Fe^2+"], ["Fe^3+"])

    # 7. 多组独立解：C + O2 + CO -> CO2
    show("多组独立解", ["C", "O2", "CO"], ["CO2"])

    # 8. 非法化学式（点号结晶水不支持；裸 +/- 电荷记号非法）
    show("非法化学式", ["CuSO4.5H2O"], ["CuSO4", "H2O"])
    show("非法电荷写法", ["Na^", "SO4^02-"], ["e^+"])

    # 9. 系数比约束：C + O2 + CO -> CO2 原本有多组独立解，
    #    约束 c(C) : c(CO) = 1 : 1 后收敛到唯一解 2C + 3O2 + 2CO -> 4CO2
    show(
        "系数比约束（多解收敛）",
        ["C", "O2", "CO"],
        ["CO2"],
        ratios=[
            {
                "left": {"side": "reactants", "index": 0},
                "right": {"side": "reactants", "index": 2},
                "ratio": [1, 1],
            }
        ],
    )

    # 10. 比值约束矛盾：H2 + O2 -> H2O 要求 2:1，强行约束 1:1 后无法配平
    show(
        "系数比约束（矛盾）",
        ["H2", "O2"],
        ["H2O"],
        ratios=[
            {
                "left": {"side": "reactants", "index": 0},
                "right": {"side": "reactants", "index": 1},
                "ratio": [1, 1],
            }
        ],
    )

    # 11. 非法比值约束（下标越界 / 同一物质 / 非正整数）：HTTP 422
    show(
        "非法比值约束",
        ["H2", "O2"],
        ["H2O"],
        ratios=[
            {
                "left": {"side": "reactants", "index": 0},
                "right": {"side": "reactants", "index": 0},
                "ratio": [1, 1],
            }
        ],
    )

    # 12. 酸性半反应自动补全：只给骨架 MnO4^- -> Mn^2+，
    #     H2O、H^+、e^- 由接口求解并放到正确的一侧
    show_acidic_half("酸性半反应（自动补全）", "MnO4^-", "Mn^2+")

    # 13. 骨架混入自动补全物种：H2O 不得出现在骨架中（HTTP 400）
    show_acidic_half("酸性半反应（骨架含 H2O）", "H2O", "Mn^2+")


if __name__ == "__main__":
    main()
