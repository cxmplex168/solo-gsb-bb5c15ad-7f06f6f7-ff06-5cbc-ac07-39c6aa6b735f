"""基于元素守恒与电荷守恒矩阵的反应配平。

把每种物质作为一列、每种元素作为一行，并额外加入一行带符号电荷：

* 反应物列填入该物质中各元素的原子数（正号）；
* 生成物列填入原子数的相反数（负号）；
* 电荷行填入各物种的带符号电荷（反应物为正、生成物为负），
  中性物种为 0，电子 ``e^-`` 为 -1。

配平系数向量 ``c`` 满足 ``A · c = 0``（原子守恒且电荷守恒）。
用 SymPy 在整数/有理数域上求精确零空间：零空间维数必须恰为 1，
且该基向量可以整体取号为全部严格正整数，再约去最大公约数得到
最小正整数解。不会为了配平而自动增删任何物种。

可选的系数比约束：``ratios`` 中每项要求两种不同物质的系数比为
``a : b``（``a``、``b`` 为正整数），转成精确整数行
``b·c_left − a·c_right = 0`` 追加到矩阵末尾，与元素、电荷守恒行
联立求解。约束可能把原本有多组独立解的反应收敛到唯一全正解；
约束矛盾或约束后不存在全正解时判为不可配平，约束后零空间维数
仍大于 1 时判为多组独立解。
"""

from __future__ import annotations

from math import gcd
from typing import NamedTuple

from sympy import Matrix, Rational

from .parser import FormulaSyntaxError, parse_species


class BalanceError(Exception):
    """配平失败基类。"""


class UnbalanceableError(BalanceError):
    """零空间无解（无法满足元素守恒与电荷守恒）。"""

    kind = "unbalanceable"


class MultipleSolutionsError(BalanceError):
    """零空间维数大于 1，存在多组独立配平。"""

    kind = "multiple_solutions"


class RatioSpec(NamedTuple):
    """一条系数比约束：``c(left) : c(right) = a : b``。

    ``left`` / ``right`` 为 ``(side, index)``：``side`` 取
    ``"reactants"`` 或 ``"products"``，``index`` 为该侧从 0 开始的下标；
    ``a``、``b`` 必须为正整数。
    """

    left: tuple[str, int]
    right: tuple[str, int]
    a: int
    b: int


def balance_reaction(
    reactants: list[str],
    products: list[str],
    ratios: list[RatioSpec] | None = None,
) -> dict:
    """配平一个化学反应（支持离子与含电子 ``e^-`` 的半反应）。

    :param ratios: 可选的系数比约束，每项作为精确矩阵行与守恒方程联立；
    :returns: 成功时返回字典，字段见模块级说明；
    :raises UnbalanceableError: 无法配平（含约束矛盾或约束后无全正解）；
    :raises MultipleSolutionsError: 约束后仍存在多组独立解；
    :raises FormulaSyntaxError: 化学式非法；
    :raises ValueError: 比值约束本身非法（侧别、越界、同一物质、非正整数）。
    """
    species = [*reactants, *products]
    n_reactants = len(reactants)
    n_products = len(products)

    parsed = [parse_species(formula) for formula in species]
    atom_counts = [counts for counts, _ in parsed]
    charges = [charge for _, charge in parsed]
    elements = sorted({el for counts in atom_counts for el in counts})

    def side_sign(index: int) -> int:
        return 1 if index < n_reactants else -1

    rows = [
        [
            Rational(counts.get(el, 0)) * side_sign(i)
            for i, counts in enumerate(atom_counts)
        ]
        for el in elements
    ]
    # 带符号电荷行：与元素行同一符号约定，保证两侧电荷守恒；
    # 全中性反应时该行为全零，不改变零空间，保持原有行为。
    rows.append(
        [Rational(charge) * side_sign(i) for i, charge in enumerate(charges)]
    )
    # 可选的系数比约束行：b·c_left − a·c_right = 0，与守恒行精确联立。
    if ratios:
        rows.extend(_ratio_rows(ratios, n_reactants, n_products))

    matrix = Matrix(rows)

    nullspace = matrix.nullspace()
    if not nullspace:
        if ratios:
            raise UnbalanceableError(
                "元素/电荷守恒与比值约束联立后无解，该反应无法配平"
            )
        raise UnbalanceableError("元素守恒方程无解，该反应无法配平")
    if len(nullspace) > 1:
        scope = "加入比值约束后" if ratios else ""
        raise MultipleSolutionsError(
            f"{scope}零空间维数为 {len(nullspace)}（大于 1），存在多组独立的配平系数"
        )

    basis = nullspace[0]
    if any(component == 0 for component in basis):
        raise UnbalanceableError(
            "唯一解中存在系数为 0 的物质，无法在保留全部输入物质的情况下配平"
        )
    if not (all(v > 0 for v in basis) or all(v < 0 for v in basis)):
        raise UnbalanceableError(
            "唯一解的系数存在正负混合，无法通过整体取号得到全部为正的系数"
        )

    coefficients = _to_coprime_positive_ints(basis)

    reactant_coeffs = coefficients[:n_reactants]
    product_coeffs = coefficients[n_reactants:]
    left_atoms = _atom_totals(atom_counts[:n_reactants], reactant_coeffs)
    right_atoms = _atom_totals(atom_counts[n_reactants:], product_coeffs)

    # 最终核验：两侧各元素原子数必须完全一致。
    if left_atoms != right_atoms:  # 理论上不会发生，留作防御
        raise UnbalanceableError("内部校验失败：两侧原子数不一致")

    left_charge = _charge_total(charges[:n_reactants], reactant_coeffs)
    right_charge = _charge_total(charges[n_reactants:], product_coeffs)

    # 最终核验：两侧带符号总电荷必须一致。
    if left_charge != right_charge:  # 理论上不会发生，留作防御
        raise UnbalanceableError("内部校验失败：两侧电荷不守恒")

    # 最终核验：每条比值约束 c_left·b = c_right·a 必须成立。
    for spec in ratios or ():  # 理论上不会发生，留作防御
        left_coeff = coefficients[
            _column_index(spec.left, n_reactants, n_products)
        ]
        right_coeff = coefficients[
            _column_index(spec.right, n_reactants, n_products)
        ]
        if left_coeff * spec.b != right_coeff * spec.a:
            raise UnbalanceableError("内部校验失败：比值约束未被满足")

    return {
        "balanced": True,
        "reactants": [
            {"formula": formula, "coefficient": coeff}
            for formula, coeff in zip(reactants, reactant_coeffs)
        ],
        "products": [
            {"formula": formula, "coefficient": coeff}
            for formula, coeff in zip(products, product_coeffs)
        ],
        "coefficients": {
            "reactants": reactant_coeffs,
            "products": product_coeffs,
        },
        "equation": _render_equation(
            reactants, reactant_coeffs, products, product_coeffs
        ),
        "atom_totals": {
            "reactants": left_atoms,
            "products": right_atoms,
        },
        "charge_totals": {
            "reactants": left_charge,
            "products": right_charge,
        },
    }


def _ratio_rows(
    ratios: list[RatioSpec], n_reactants: int, n_products: int
) -> list[list[Rational]]:
    """把每条比值约束 ``c_left : c_right = a : b`` 转成精确整数行
    ``b·c_left − a·c_right = 0``（其余列为 0）。"""
    n_species = n_reactants + n_products
    rows: list[list[Rational]] = []
    for spec in ratios:
        if spec.a <= 0 or spec.b <= 0:
            raise ValueError("比值约束的 a、b 必须为正整数")
        left_col = _column_index(spec.left, n_reactants, n_products)
        right_col = _column_index(spec.right, n_reactants, n_products)
        if left_col == right_col:
            raise ValueError("比值约束的 left 与 right 必须指向两种不同物质")
        row = [Rational(0)] * n_species
        row[left_col] = Rational(spec.b)
        row[right_col] = Rational(-spec.a)
        rows.append(row)
    return rows


def _column_index(
    ref: tuple[str, int], n_reactants: int, n_products: int
) -> int:
    """把 ``(side, index)`` 引用解析为矩阵列号（反应物列在前）。"""
    side, index = ref
    if side == "reactants":
        base, size = 0, n_reactants
    elif side == "products":
        base, size = n_reactants, n_products
    else:
        raise ValueError(f"未知的侧别 '{side}'，应为 'reactants' 或 'products'")
    if not 0 <= index < size:
        raise ValueError(
            f"{side} 侧下标 {index} 越界（该侧共 {size} 种物质）"
        )
    return base + index


def _to_coprime_positive_ints(basis: Matrix) -> list[int]:
    """把有理零空间基向量取号为正并约成 gcd = 1 的最小正整数。"""
    rationals = [Rational(v) for v in basis]
    denominators = [rational.q for rational in rationals]

    common_denominator = 1
    for denominator in denominators:
        common_denominator = common_denominator * denominator // gcd(
            common_denominator, denominator
        )

    values = [int(rational * common_denominator) for rational in rationals]
    if all(v < 0 for v in values):
        values = [-v for v in values]

    common_divisor = 0
    for value in values:
        common_divisor = gcd(common_divisor, abs(value))
    if common_divisor > 1:
        values = [v // common_divisor for v in values]
    return values


def _atom_totals(
    parsed_side: list[dict[str, int]], coeffs: list[int]
) -> dict[str, int]:
    totals: dict[str, int] = {}
    for counts, coeff in zip(parsed_side, coeffs):
        for element, amount in counts.items():
            totals[element] = totals.get(element, 0) + amount * coeff
    return dict(sorted(totals.items()))


def _charge_total(charges: list[int], coeffs: list[int]) -> int:
    """一侧的带符号总电荷：Σ 系数 × 物种电荷。"""
    return sum(charge * coeff for charge, coeff in zip(charges, coeffs))


def _render_equation(
    reactants: list[str],
    reactant_coeffs: list[int],
    products: list[str],
    product_coeffs: list[int],
) -> str:
    def side(formulas: list[str], coeffs: list[int]) -> str:
        parts = []
        for formula, coeff in zip(formulas, coeffs):
            parts.append(formula if coeff == 1 else f"{coeff}{formula}")
        return " + ".join(parts)

    return f"{side(reactants, reactant_coeffs)} -> {side(products, product_coeffs)}"
