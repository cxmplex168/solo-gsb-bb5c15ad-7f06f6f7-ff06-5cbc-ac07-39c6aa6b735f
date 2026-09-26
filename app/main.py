"""FastAPI 入口。"""

from __future__ import annotations

import os

from fastapi import FastAPI
from fastapi.responses import JSONResponse

from .balancer import (
    ACIDIC_AUXILIARY_SPECIES,
    MultipleSolutionsError,
    NoElectronTransferError,
    RatioSpec,
    UnbalanceableError,
    balance_acidic_half_reaction,
    balance_reaction,
)
from .parser import FormulaSyntaxError, parse_species
from .schemas import AcidicHalfRequest, BalanceRequest

app = FastAPI(
    title="化学反应配平 API",
    description=(
        "输入反应物与生成物的化学式数组，基于元素守恒与电荷守恒的精确整数矩阵配平。"
        "支持元素符号、正整数下标、最多四层圆括号、末尾 ^+/^-/^n+/^n- 电荷后缀"
        "与电子 e^-（可配平离子反应与半反应）；不支持点号结晶水。"
        "可选 ratios 系数比约束与守恒方程联立求解，可把多解反应收敛到唯一解。"
        "/balance/acidic-half 接收仅含一种反应物和一种生成物的半反应骨架，"
        "在酸性条件下自动补全 H2O、H^+、e^- 并配平。"
    ),
    version="1.3.0",
)


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.post(
    "/balance",
    summary="配平化学反应",
    responses={
        200: {
            "description": (
                "配平结果或明确的失败原因。注意：无法配平 / 多组独立解 / "
                "非法化学式以 HTTP 400 或 200 返回，通过 balanced 与 reason 字段区分；"
                "ratios 约束的非法侧别、下标越界、同一物质或非正整数以 HTTP 422 返回。"
            )
        }
    },
)
def balance(request: BalanceRequest) -> JSONResponse:
    # 先统一校验两侧所有化学式，一次性返回全部非法位置。
    formula_errors: list[dict] = []
    for side_name, formulas in (
        ("reactants", request.reactants),
        ("products", request.products),
    ):
        for index, formula in enumerate(formulas):
            try:
                parse_species(formula)
            except FormulaSyntaxError as exc:
                formula_errors.append(
                    {
                        "side": side_name,
                        "index": index,
                        "formula": formula,
                        "position": exc.position,
                        "message": exc.message,
                    }
                )

    if formula_errors:
        return JSONResponse(
            status_code=400,
            content={
                "balanced": False,
                "reason": "invalid_formula",
                "message": "存在非法化学式，错误位置见 errors",
                "errors": formula_errors,
            },
        )

    # 把请求中的比值约束转成求解器使用的 (side, index) 引用形式；
    # 侧别、越界、同一物质、非正整数已由请求模型校验（HTTP 422）。
    ratio_specs = [
        RatioSpec(
            left=(constraint.left.side, constraint.left.index),
            right=(constraint.right.side, constraint.right.index),
            a=constraint.ratio[0],
            b=constraint.ratio[1],
        )
        for constraint in request.ratios or ()
    ]

    try:
        result = balance_reaction(
            request.reactants, request.products, ratios=ratio_specs
        )
    except UnbalanceableError as exc:
        return JSONResponse(
            status_code=200,
            content={
                "balanced": False,
                "reason": "unbalanceable",
                "message": str(exc),
            },
        )
    except MultipleSolutionsError as exc:
        return JSONResponse(
            status_code=200,
            content={
                "balanced": False,
                "reason": "multiple_solutions",
                "message": str(exc),
            },
        )

    return JSONResponse(status_code=200, content=result)


def _species_key(counts: dict[str, int], charge: int) -> tuple:
    """物种解析结果的规范化键，用于识别语义上相同的化学式。"""
    return tuple(sorted(counts.items())), charge


# 酸性半反应接口中由系统自动补全、禁止出现在骨架里的物种
# （按解析结果匹配，HOH、H1^+ 等等价写法同样会被拒绝）
_FORBIDDEN_SKELETON_SPECIES = {
    _species_key(*parse_species(formula)): formula
    for formula in ACIDIC_AUXILIARY_SPECIES
}


@app.post(
    "/balance/acidic-half",
    summary="酸性条件下自动补全并配平半反应",
    responses={
        200: {
            "description": (
                "配平结果或明确的失败原因。无法配平 / 多组独立解 / 电子系数为 0 "
                "以 HTTP 200 返回，通过 balanced 与 reason 字段区分；"
                "非法化学式或骨架含 H2O、H^+、e^- 以 HTTP 400 返回；"
                "任一侧化学式数量不为 1 以 HTTP 422 返回。"
            )
        }
    },
)
def balance_acidic_half(request: AcidicHalfRequest) -> JSONResponse:
    # 先统一校验两侧化学式，并检查是否混入自动补全物种（H2O、H^+、e^-）。
    formula_errors: list[dict] = []
    forbidden_species: list[dict] = []
    for side_name, formulas in (
        ("reactants", request.reactants),
        ("products", request.products),
    ):
        for index, formula in enumerate(formulas):
            try:
                counts, charge = parse_species(formula)
            except FormulaSyntaxError as exc:
                formula_errors.append(
                    {
                        "side": side_name,
                        "index": index,
                        "formula": formula,
                        "position": exc.position,
                        "message": exc.message,
                    }
                )
                continue
            canonical = _FORBIDDEN_SKELETON_SPECIES.get(
                _species_key(counts, charge)
            )
            if canonical is not None:
                forbidden_species.append(
                    {
                        "side": side_name,
                        "index": index,
                        "formula": formula,
                        "message": (
                            f"{canonical} 由接口在酸性条件下自动补全，"
                            "不得出现在半反应骨架中"
                        ),
                    }
                )

    if formula_errors:
        return JSONResponse(
            status_code=400,
            content={
                "balanced": False,
                "reason": "invalid_formula",
                "message": "存在非法化学式，错误位置见 errors",
                "errors": formula_errors,
            },
        )
    if forbidden_species:
        return JSONResponse(
            status_code=400,
            content={
                "balanced": False,
                "reason": "forbidden_species",
                "message": (
                    "骨架中不得包含 H2O、H^+、e^-（由接口自动补全），"
                    "详见 errors"
                ),
                "errors": forbidden_species,
            },
        )

    try:
        result = balance_acidic_half_reaction(
            request.reactants[0], request.products[0]
        )
    except UnbalanceableError as exc:
        return JSONResponse(
            status_code=200,
            content={
                "balanced": False,
                "reason": "unbalanceable",
                "message": str(exc),
            },
        )
    except MultipleSolutionsError as exc:
        return JSONResponse(
            status_code=200,
            content={
                "balanced": False,
                "reason": "multiple_solutions",
                "message": str(exc),
            },
        )
    except NoElectronTransferError as exc:
        return JSONResponse(
            status_code=200,
            content={
                "balanced": False,
                "reason": "no_electron_transfer",
                "message": str(exc),
            },
        )

    return JSONResponse(status_code=200, content=result)


if __name__ == "__main__":
    import uvicorn

    host = os.environ.get("APP_HOST", "0.0.0.0")
    port = int(os.environ.get("APP_PORT", "8000"))
    uvicorn.run("app.main:app", host=host, port=port)
