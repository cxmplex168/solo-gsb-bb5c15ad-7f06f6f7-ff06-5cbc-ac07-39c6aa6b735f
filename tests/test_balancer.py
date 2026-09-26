import pytest

from app.balancer import (
    MultipleSolutionsError,
    NoElectronTransferError,
    RatioSpec,
    UnbalanceableError,
    balance_acidic_half_reaction,
    balance_reaction,
)


def test_simple_balance():
    result = balance_reaction(["Fe", "O2"], ["Fe3O4"])
    assert result["coefficients"] == {"reactants": [3, 2], "products": [1]}
    assert result["equation"] == "3Fe + 2O2 -> Fe3O4"
    assert result["atom_totals"]["reactants"] == {"Fe": 3, "O": 4}
    assert result["atom_totals"]["products"] == {"Fe": 3, "O": 4}


def test_order_is_preserved():
    result = balance_reaction(["O2", "Fe"], ["Fe3O4"])
    assert result["coefficients"] == {"reactants": [2, 3], "products": [1]}
    assert result["reactants"][0]["formula"] == "O2"
    assert result["reactants"][1]["formula"] == "Fe"


def test_water_electrolysis():
    result = balance_reaction(["H2O"], ["H2", "O2"])
    assert result["coefficients"] == {"reactants": [2], "products": [2, 1]}
    totals = result["atom_totals"]
    assert totals["reactants"] == totals["products"] == {"H": 4, "O": 2}


def test_coefficients_are_minimal_gcd_one():
    # KMnO4 + HCl -> KCl + MnCl2 + Cl2 + H2O 的经典配平
    result = balance_reaction(
        ["KMnO4", "HCl"],
        ["KCl", "MnCl2", "Cl2", "H2O"],
    )
    coeffs = (
        result["coefficients"]["reactants"]
        + result["coefficients"]["products"]
    )
    assert coeffs == [2, 16, 2, 2, 5, 8]
    from math import gcd
    from functools import reduce

    assert reduce(gcd, coeffs) == 1


def test_nested_parentheses_balance():
    # 黄血盐 + 高锰酸钾酸性体系（示例中使用的反应）
    result = balance_reaction(
        ["K4(Fe(CN)6)", "KMnO4", "H2SO4"],
        ["KHSO4", "Fe2(SO4)3", "MnSO4", "HNO3", "CO2", "H2O"],
    )
    left = result["atom_totals"]["reactants"]
    right = result["atom_totals"]["products"]
    assert left == right
    assert all(c > 0 for c in result["coefficients"]["reactants"])
    assert all(c > 0 for c in result["coefficients"]["products"])
    assert result["coefficients"]["reactants"][0] == 10
    assert result["coefficients"]["reactants"][1] == 122


def test_unbalanceable_full_rank_no_solution():
    # 产物多了一个反应物中不存在的元素，且约束满秩：零空间为空
    with pytest.raises(UnbalanceableError):
        balance_reaction(["H2", "O2"], ["HFeO2"])


def test_unbalanceable_zero_coefficient_distinguished():
    # H2 + O2 -> H2O + Fe：Fe 列被迫为 0，唯一解含零系数
    with pytest.raises(UnbalanceableError):
        balance_reaction(["H2", "O2"], ["H2O", "Fe"])


def test_multiple_solutions():
    # C + O2 + CO -> CO2：CO + 1/2 O2 -> CO2 与 C + O2 -> CO2 两个独立反应
    with pytest.raises(MultipleSolutionsError):
        balance_reaction(["C", "O2", "CO"], ["CO2"])


def test_zero_coefficient_solution_is_unbalanceable():
    # 体系中混入一个与反应无关的物质时，唯一解里它的系数可能是 0
    with pytest.raises(UnbalanceableError):
        balance_reaction(["H2", "O2", "N2"], ["H2O"])


def test_single_species_reaction():
    # 同素异形体：O3 -> O2
    result = balance_reaction(["O3"], ["O2"])
    assert result["coefficients"] == {"reactants": [2], "products": [3]}


def test_repeated_formula_entries():
    # 同一物质出现两次：H2 + H2 -> H2O ... 无法配平（O 无来源）
    with pytest.raises(UnbalanceableError):
        balance_reaction(["H2", "H2"], ["H2O"])


def test_neutral_reaction_has_zero_charge_totals():
    result = balance_reaction(["Fe", "O2"], ["Fe3O4"])
    assert result["charge_totals"] == {"reactants": 0, "products": 0}


def test_half_reaction_with_electron():
    # 氧化半反应：Fe^2+ -> Fe^3+ + e^-
    result = balance_reaction(["Fe^2+"], ["Fe^3+", "e^-"])
    assert result["coefficients"] == {"reactants": [1], "products": [1, 1]}
    assert result["charge_totals"] == {"reactants": 2, "products": 2}
    assert result["atom_totals"]["reactants"] == {"Fe": 1}
    assert result["atom_totals"]["products"] == {"Fe": 1}


def test_acidic_permanganate_half_reaction():
    # 酸性条件下的还原半反应：MnO4^- + H^+ + e^- -> Mn^2+ + H2O
    result = balance_reaction(
        ["MnO4^-", "H^+", "e^-"],
        ["Mn^2+", "H2O"],
    )
    assert result["coefficients"] == {"reactants": [1, 8, 5], "products": [1, 4]}
    assert result["charge_totals"] == {"reactants": 2, "products": 2}
    totals = result["atom_totals"]
    assert totals["reactants"] == totals["products"] == {"H": 8, "Mn": 1, "O": 4}
    assert result["equation"] == "MnO4^- + 8H^+ + 5e^- -> Mn^2+ + 4H2O"


def test_ionic_reaction_without_electron():
    # 离子方程式：Zn + Cu^2+ -> Zn^2+ + Cu
    result = balance_reaction(["Zn", "Cu^2+"], ["Zn^2+", "Cu"])
    assert result["coefficients"] == {"reactants": [1, 1], "products": [1, 1]}
    assert result["charge_totals"] == {"reactants": 2, "products": 2}


def test_charge_not_conserved_is_unbalanceable():
    # 原子守恒但电荷不守恒：Fe^2+ -> Fe^3+ 缺少电子
    with pytest.raises(UnbalanceableError):
        balance_reaction(["Fe^2+"], ["Fe^3+"])


def test_electron_only_on_one_side_is_unbalanceable():
    # 电子被迫系数为 0：H2O -> H2 + O2 + e^-
    with pytest.raises(UnbalanceableError):
        balance_reaction(["H2O"], ["H2", "O2", "e^-"])


def test_ratio_constraint_converges_multiple_to_unique():
    # C + O2 + CO -> CO2 原本有两组独立解；约束 c(C) : c(CO) = 1 : 1 后收敛
    result = balance_reaction(
        ["C", "O2", "CO"],
        ["CO2"],
        ratios=[RatioSpec(left=("reactants", 0), right=("reactants", 2), a=1, b=1)],
    )
    assert result["coefficients"] == {"reactants": [2, 3, 2], "products": [4]}
    assert result["equation"] == "2C + 3O2 + 2CO -> 4CO2"
    totals = result["atom_totals"]
    assert totals["reactants"] == totals["products"] == {"C": 4, "O": 8}


def test_ratio_constraint_across_sides():
    # 约束反应物 H2 与生成物 H2O 的系数比为 1:1
    result = balance_reaction(
        ["H2", "O2"],
        ["H2O"],
        ratios=[RatioSpec(("reactants", 0), ("products", 0), 1, 1)],
    )
    assert result["coefficients"] == {"reactants": [2, 1], "products": [2]}


def test_ratio_constraint_non_coprime_ratio_is_equivalent():
    # 2:2 与 1:1 等效，结果仍约成 gcd = 1 的最小整数
    result = balance_reaction(
        ["C", "O2", "CO"],
        ["CO2"],
        ratios=[RatioSpec(("reactants", 0), ("reactants", 2), 2, 2)],
    )
    coeffs = (
        result["coefficients"]["reactants"] + result["coefficients"]["products"]
    )
    assert coeffs == [2, 3, 2, 4]
    from functools import reduce
    from math import gcd

    assert reduce(gcd, coeffs) == 1


def test_ratio_constraint_contradiction_is_unbalanceable():
    # H2 + O2 -> H2O 要求 2:1，强行约束 c(H2) : c(O2) = 1 : 1 后联立无解
    with pytest.raises(UnbalanceableError):
        balance_reaction(
            ["H2", "O2"],
            ["H2O"],
            ratios=[RatioSpec(("reactants", 0), ("reactants", 1), 1, 1)],
        )


def test_ratio_constraint_zero_coefficient_is_unbalanceable():
    # 约束收敛到唯一解，但与反应无关的 N2 系数被迫为 0
    with pytest.raises(UnbalanceableError):
        balance_reaction(
            ["H2", "O2", "N2"],
            ["H2O"],
            ratios=[RatioSpec(("reactants", 0), ("reactants", 1), 2, 1)],
        )


def test_ratio_constraint_still_multiple_solutions():
    # 零空间原为 3 维，一条约束只降到 2 维，仍为多组独立解
    with pytest.raises(MultipleSolutionsError):
        balance_reaction(
            ["H2", "O2", "C", "CO"],
            ["H2O", "CO2"],
            ratios=[RatioSpec(("reactants", 0), ("reactants", 1), 2, 1)],
        )


def test_ratio_omitted_or_empty_keeps_original_behavior():
    plain = balance_reaction(["Fe", "O2"], ["Fe3O4"])
    empty = balance_reaction(["Fe", "O2"], ["Fe3O4"], ratios=[])
    expected = {"reactants": [3, 2], "products": [1]}
    assert plain["coefficients"] == empty["coefficients"] == expected


def test_ratio_invalid_spec_raises_value_error():
    # 越界
    with pytest.raises(ValueError):
        balance_reaction(
            ["H2", "O2"],
            ["H2O"],
            ratios=[RatioSpec(("reactants", 0), ("reactants", 2), 1, 1)],
        )
    # 同一物质
    with pytest.raises(ValueError):
        balance_reaction(
            ["H2", "O2"],
            ["H2O"],
            ratios=[RatioSpec(("products", 0), ("products", 0), 1, 1)],
        )
    # 非正整数
    with pytest.raises(ValueError):
        balance_reaction(
            ["H2", "O2"],
            ["H2O"],
            ratios=[RatioSpec(("reactants", 0), ("reactants", 1), 0, 1)],
        )
    # 未知侧别
    with pytest.raises(ValueError):
        balance_reaction(
            ["H2", "O2"],
            ["H2O"],
            ratios=[RatioSpec(("left", 0), ("reactants", 1), 1, 1)],
        )


def test_acidic_half_permanganate_reduction():
    # 经典酸性还原半反应：MnO4^- + 8H^+ + 5e^- -> Mn^2+ + 4H2O
    result = balance_acidic_half_reaction("MnO4^-", "Mn^2+")
    assert result["coefficients"] == {"reactants": [1, 8, 5], "products": [1, 4]}
    assert result["equation"] == "MnO4^- + 8H^+ + 5e^- -> Mn^2+ + 4H2O"
    assert result["charge_totals"] == {"reactants": 2, "products": 2}
    totals = result["atom_totals"]
    assert totals["reactants"] == totals["products"] == {"H": 8, "Mn": 1, "O": 4}
    # 骨架物质排在各自侧最前，自动补全项按 H2O、H^+、e^- 顺序追加
    assert [s["formula"] for s in result["reactants"]] == ["MnO4^-", "H^+", "e^-"]
    assert [s["formula"] for s in result["products"]] == ["Mn^2+", "H2O"]


def test_acidic_half_zero_terms_dropped():
    # Fe^2+ -> Fe^3+ + e^-：H2O 与 H^+ 系数为 0，直接删去
    result = balance_acidic_half_reaction("Fe^2+", "Fe^3+")
    assert result["coefficients"] == {"reactants": [1], "products": [1, 1]}
    assert result["equation"] == "Fe^2+ -> Fe^3+ + e^-"
    assert result["charge_totals"] == {"reactants": 2, "products": 2}


def test_acidic_half_water_on_reactant_side():
    # H2O 解为负值移到反应物侧，H^+ 与 e^- 为正值放到生成物侧
    result = balance_acidic_half_reaction("SO3^2-", "SO4^2-")
    assert result["coefficients"] == {"reactants": [1, 1], "products": [1, 2, 2]}
    assert result["equation"] == "SO3^2- + H2O -> SO4^2- + 2H^+ + 2e^-"
    assert result["charge_totals"] == {"reactants": -2, "products": -2}
    totals = result["atom_totals"]
    assert totals["reactants"] == totals["products"] == {"H": 2, "O": 4, "S": 1}


def test_acidic_half_water_on_product_side():
    # MnO2 + 4H^+ + 2e^- -> Mn^2+ + 2H2O
    result = balance_acidic_half_reaction("MnO2", "Mn^2+")
    assert result["coefficients"] == {"reactants": [1, 4, 2], "products": [1, 2]}
    assert result["equation"] == "MnO2 + 4H^+ + 2e^- -> Mn^2+ + 2H2O"


def test_acidic_half_dichromate_coprime_coefficients():
    # Cr2O7^2- + 14H^+ + 6e^- -> 2Cr^3+ + 7H2O，全体系数互质
    result = balance_acidic_half_reaction("Cr2O7^2-", "Cr^3+")
    coeffs = (
        result["coefficients"]["reactants"] + result["coefficients"]["products"]
    )
    assert coeffs == [1, 14, 6, 2, 7]
    from functools import reduce
    from math import gcd

    assert reduce(gcd, coeffs) == 1


def test_acidic_half_electron_on_product_side():
    # 氧化半反应：C2O4^2- -> 2CO2 + 2e^-（H2O、H^+ 系数为 0 被删去）
    result = balance_acidic_half_reaction("C2O4^2-", "CO2")
    assert result["coefficients"] == {"reactants": [1], "products": [2, 2]}
    assert result["equation"] == "C2O4^2- -> 2CO2 + 2e^-"
    assert result["charge_totals"] == {"reactants": -2, "products": -2}


def test_acidic_half_no_electron_transfer():
    # H2CO3 -> HCO3^- + H^+ 是酸碱解离，电子系数为 0
    with pytest.raises(NoElectronTransferError):
        balance_acidic_half_reaction("H2CO3", "HCO3^-")


def test_acidic_half_unbalanceable_no_solution():
    # K 元素无法守恒：KMnO4 -> K2MnO4 零空间为空
    with pytest.raises(UnbalanceableError):
        balance_acidic_half_reaction("KMnO4", "K2MnO4")


def test_acidic_half_unbalanceable_zero_skeleton_coefficient():
    # Fe 与 O2 无共同元素，唯一解中骨架反应物系数为 0
    with pytest.raises(UnbalanceableError):
        balance_acidic_half_reaction("Fe", "O2")


def test_acidic_half_multiple_solutions():
    # 骨架只含 H、O 元素时 H2O / H^+ 无法唯一确定：零空间维数为 2
    with pytest.raises(MultipleSolutionsError):
        balance_acidic_half_reaction("H2O2", "O2")
