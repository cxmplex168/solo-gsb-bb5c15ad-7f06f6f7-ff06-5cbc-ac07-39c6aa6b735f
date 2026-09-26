"""请求/响应模型。"""

from __future__ import annotations

from typing import Annotated, Literal

from pydantic import BaseModel, Field, model_validator

FormulaList = Annotated[
    list[str],
    Field(
        min_length=1,
        max_length=50,
        description="化学式数组，保持数组顺序作为返回系数顺序",
    ),
]

SingleFormula = Annotated[
    list[str],
    Field(
        min_length=1,
        max_length=1,
        description="恰好一个化学式的数组（半反应骨架的一侧）",
    ),
]

Side = Literal["reactants", "products"]

PositiveInt = Annotated[int, Field(gt=0, description="正整数")]


class SpeciesRef(BaseModel):
    """对某一侧某一种物质的引用（侧别 + 该侧 0 起始下标）。"""

    side: Side = Field(description="物质所在侧：reactants 或 products")
    index: int = Field(ge=0, description="该侧数组中的下标（从 0 开始）")


class RatioConstraint(BaseModel):
    """系数比约束：c(left) : c(right) = ratio[0] : ratio[1]。"""

    left: SpeciesRef
    right: SpeciesRef
    ratio: tuple[PositiveInt, PositiveInt] = Field(
        description="[a, b] 两个正整数，要求 left 与 right 所指物质的系数比为 a:b"
    )


class BalanceRequest(BaseModel):
    reactants: FormulaList
    products: FormulaList
    ratios: Annotated[
        list[RatioConstraint],
        Field(
            max_length=100,
            description=(
                "可选的系数比约束数组；每项要求两种不同物质的系数比为 a:b，"
                "与元素、电荷守恒联立求解。省略时保持原有配平行为"
            ),
        ),
    ] | None = None

    @model_validator(mode="after")
    def _check_ratio_refs(self) -> "BalanceRequest":
        """校验比值约束的下标越界与同一物质引用（返回 HTTP 422）。"""
        if not self.ratios:
            return self
        sizes = {"reactants": len(self.reactants), "products": len(self.products)}
        for pos, constraint in enumerate(self.ratios):
            for label, ref in (
                ("left", constraint.left),
                ("right", constraint.right),
            ):
                size = sizes[ref.side]
                if ref.index >= size:
                    raise ValueError(
                        f"ratios[{pos}].{label}: {ref.side} 侧下标 {ref.index} 越界，"
                        f"该侧共 {size} 种物质（合法下标 0 到 {size - 1}）"
                    )
            if (
                constraint.left.side == constraint.right.side
                and constraint.left.index == constraint.right.index
            ):
                raise ValueError(
                    f"ratios[{pos}]: left 与 right 指向同一物质 "
                    f"{constraint.left.side}[{constraint.left.index}]，"
                    "比值约束必须作用于两种不同物质"
                )
        return self


class AcidicHalfBalanceRequest(BaseModel):
    """酸性半反应配平请求：每侧恰好一个骨架化学式。

    骨架不得包含 H2O、H^+、e^-（作为 HTTP 400 业务校验，由接口处理）。
    """

    reactants: SingleFormula
    products: SingleFormula


class SpeciesCoefficient(BaseModel):
    formula: str
    coefficient: int


class CoefficientPair(BaseModel):
    reactants: list[int]
    products: list[int]


class AtomTotals(BaseModel):
    reactants: dict[str, int]
    products: dict[str, int]


class ChargeTotals(BaseModel):
    reactants: int = Field(description="反应物侧带符号总电荷")
    products: int = Field(description="生成物侧带符号总电荷")


class BalanceSuccessResponse(BaseModel):
    balanced: Literal[True] = True
    reactants: list[SpeciesCoefficient]
    products: list[SpeciesCoefficient]
    coefficients: CoefficientPair
    equation: str
    atom_totals: AtomTotals
    charge_totals: ChargeTotals


class FormulaErrorDetail(BaseModel):
    side: Literal["reactants", "products"]
    index: int = Field(description="该侧数组中的下标（从 0 开始）")
    formula: str
    position: int = Field(description="化学式中的非法字符位置（从 0 开始）")
    message: str


class ForbiddenSpeciesDetail(BaseModel):
    side: Literal["reactants", "products"]
    index: int = Field(description="该侧数组中的下标（从 0 开始）")
    formula: str


class BalanceFailureResponse(BaseModel):
    balanced: Literal[False] = False
    reason: Literal[
        "invalid_formula",
        "forbidden_species",
        "unbalanceable",
        "multiple_solutions",
        "electron_coefficient_zero",
    ]
    message: str
    errors: list[FormulaErrorDetail | ForbiddenSpeciesDetail] | None = None
