"""酸性半反应配平 balance_acidic_half_reaction 的单元测试。"""

import pytest

from app.balancer import (
    ElectronCoefficientZeroError,
    MultipleSolutionsError,
    UnbalanceableError,
    balance_acidic_half_reaction,
)


def test_permanganate_reduction():
    # MnO4^- + 8H^+ + 5e^- -> Mn^2+ + 4H2O
    result = balance_acidic_half_reaction("MnO4^-", "Mn^2+")
    assert result["coefficients"] == {"reactants": [1, 8, 5], "products": [1, 4]}
    assert result["reactants"] == [
        {"formula": "MnO4^-", "coefficient": 1},
        {"formula": "H^+", "coefficient": 8},
        {"formula": "e^-", "coefficient": 5},
    ]
    assert result["products"] == [
        {"formula": "Mn^2+", "coefficient": 1},
        {"formula": "H2O", "coefficient": 4},
    ]
    assert result["equation"] == "MnO4^- + 8H^+ + 5e^- -> Mn^2+ + 4H2O"
    assert result["charge_totals"] == {"reactants": 2, "products": 2}
    assert result["atom_totals"] == {
        "reactants": {"H": 8, "Mn": 1, "O": 4},
        "products": {"H": 8, "Mn": 1, "O": 4},
    }


def test_dichromate_reduction():
    # Cr2O7^2- + 14H^+ + 6e^- -> 2Cr^3+ + 7H2O
    result = balance_acidic_half_reaction("Cr2O7^2-", "Cr^3+")
    assert result["coefficients"] == {"reactants": [1, 14, 6], "products": [2, 7]}
    assert result["equation"] == (
        "Cr2O7^2- + 14H^+ + 6e^- -> 2Cr^3+ + 7H2O"
    )
    assert result["charge_totals"] == {"reactants": 6, "products": 6}


def test_nitrate_to_no_reduction():
    result = balance_acidic_half_reaction("NO3^-", "NO")
    assert result["equation"] == "NO3^- + 4H^+ + 3e^- -> NO + 2H2O"
    assert result["coefficients"] == {"reactants": [1, 4, 3], "products": [1, 2]}


def test_simple_oxidation_electron_on_product_side():
    # Fe^2+ -> Fe^3+ + e^-
    result = balance_acidic_half_reaction("Fe^2+", "Fe^3+")
    assert result["coefficients"] == {"reactants": [1], "products": [1, 1]}
    assert result["products"][-1] == {"formula": "e^-", "coefficient": 1}
    assert "H2O" not in result["equation"]
    assert "H^+" not in result["equation"]
    assert result["equation"] == "Fe^2+ -> Fe^3+ + e^-"


def test_simple_reduction_electron_on_reactant_side():
    # Fe^3+ + e^- -> Fe^2+
    result = balance_acidic_half_reaction("Fe^3+", "Fe^2+")
    assert result["coefficients"] == {"reactants": [1, 1], "products": [1]}
    assert result["reactants"][-1] == {"formula": "e^-", "coefficient": 1}
    assert result["equation"] == "Fe^3+ + e^- -> Fe^2+"


def test_water_moves_to_reactant_side_in_oxidation():
    # NH4^+ + 3H2O -> NO3^- + 10H^+ + 8e^-（水为负未知量，移到反应物侧）
    result = balance_acidic_half_reaction("NH4^+", "NO3^-")
    assert result["coefficients"] == {"reactants": [1, 3], "products": [1, 10, 8]}
    assert result["equation"] == "NH4^+ + 3H2O -> NO3^- + 10H^+ + 8e^-"
    assert result["charge_totals"] == {"reactants": 1, "products": 1}
    assert result["atom_totals"] == {
        "reactants": {"H": 10, "N": 1, "O": 3},
        "products": {"H": 10, "N": 1, "O": 3},
    }


def test_neutral_reactant_metal_to_ion():
    # Fe -> Fe^3+ + 3e^-
    result = balance_acidic_half_reaction("Fe", "Fe^3+")
    assert result["equation"] == "Fe -> Fe^3+ + 3e^-"
    assert result["charge_totals"] == {"reactants": 0, "products": 0}


def test_oxide_skeleton_coefficients_minimal():
    # Fe2O3 + 6H^+ + 2e^- -> 2Fe^2+ + 3H2O
    result = balance_acidic_half_reaction("Fe2O3", "Fe^2+")
    assert result["coefficients"] == {"reactants": [1, 6, 2], "products": [2, 3]}
    from functools import reduce
    from math import gcd

    all_coeffs = (
        result["coefficients"]["reactants"] + result["coefficients"]["products"]
    )
    assert reduce(gcd, all_coeffs) == 1


def test_peroxodisulfate_reduction_no_hydrogen_terms():
    # S2O8^2- + 2e^- -> 2SO4^2-：H2O 与 H^+ 系数为零被删去
    result = balance_acidic_half_reaction("S2O8^2-", "SO4^2-")
    assert result["coefficients"] == {"reactants": [1, 2], "products": [2]}
    assert result["equation"] == "S2O8^2- + 2e^- -> 2SO4^2-"
    assert "H2O" not in result["equation"]
    assert "H^+" not in result["equation"]


def test_electron_coefficient_zero_is_distinct_failure():
    # P4O10 + 12H^+ -> 4H3PO4：纯水解反应，电子系数为 0
    with pytest.raises(ElectronCoefficientZeroError):
        balance_acidic_half_reaction("P4O10", "H3PO4")
    # 相同物质 Fe -> Fe：唯一解电子为 0
    with pytest.raises(ElectronCoefficientZeroError):
        balance_acidic_half_reaction("Fe", "Fe")


def test_multiple_solutions_non_unique_skeleton():
    # H2O2 -> O2：H/O 仅经 H2O/H^+ 耦合，存在多组独立半反应方向
    with pytest.raises(MultipleSolutionsError):
        balance_acidic_half_reaction("H2O2", "O2")
    with pytest.raises(MultipleSolutionsError):
        balance_acidic_half_reaction("O2", "H2O2")
    with pytest.raises(MultipleSolutionsError):
        balance_acidic_half_reaction("H2", "O2")
    # O3 -> O2：两种骨架系数都自由，同样不存在唯一全正骨架解
    with pytest.raises(MultipleSolutionsError):
        balance_acidic_half_reaction("O3", "O2")


def test_unbalanceable_skeleton():
    # 骨架元素构成无法与 H2O/H^+/e^- 联立守恒：Na 元素无对应去向
    with pytest.raises(UnbalanceableError):
        balance_acidic_half_reaction("NaCl", "Cl^-")


@pytest.mark.parametrize("forbidden", ["H2O", "H^+", "e^-"])
def test_forbidden_auxiliary_species_in_skeleton(forbidden):
    with pytest.raises(ValueError):
        balance_acidic_half_reaction(forbidden, "Fe^2+")
    with pytest.raises(ValueError):
        balance_acidic_half_reaction("Fe^3+", forbidden)


def test_invalid_formula_propagates_syntax_error():
    from app.parser import FormulaSyntaxError

    with pytest.raises(FormulaSyntaxError):
        balance_acidic_half_reaction("Na+", "Na")
