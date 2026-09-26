"""FastAPI 入口。"""

from __future__ import annotations

import os

from fastapi import FastAPI
from fastapi.responses import JSONResponse

from .balancer import (
    ACIDIC_AUX_SPECIES,
    ElectronCoefficientZeroError,
    MultipleSolutionsError,
    RatioSpec,
    UnbalanceableError,
    balance_acidic_half_reaction,
    balance_reaction,
)
from .parser import FormulaSyntaxError, parse_species
from .schemas import AcidicHalfBalanceRequest, BalanceRequest

app = FastAPI(
    title="化学反应配平 API",
    description=(
        "输入反应物与生成物的化学式数组，基于元素守恒与电荷守恒的精确整数矩阵配平。"
        "支持元素符号、正整数下标、最多四层圆括号、末尾 ^+/^-/^n+/^n- 电荷后缀"
        "与电子 e^-（可配平离子反应与半反应）；不支持点号结晶水。"
        "可选 ratios 系数比约束与守恒方程联立求解，可把多解反应收敛到唯一解。"
        "POST /balance/acidic-half 接收仅含一种反应物与一种生成物的酸性半反应骨架，"
        "自动精确补入 H2O、H^+、e^- 并配平。"
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


@app.post(
    "/balance/acidic-half",
    summary="配平酸性条件下的半反应（骨架自动补 H2O/H^+/e^-）",
    responses={
        200: {
            "description": (
                "配平结果或明确的失败原因。无法配平 / 多组独立解 / "
                "电子系数为零 / 非法化学式 / 骨架包含禁止物种均以 "
                "HTTP 200 或 400 返回，通过 balanced 与 reason 字段区分；"
                "每侧物质数不是恰好一个以 HTTP 422 返回。"
            )
        }
    },
)
def balance_acidic_half(request: AcidicHalfBalanceRequest) -> JSONResponse:
    reactant = request.reactants[0]
    product = request.products[0]

    # 先校验两个骨架化学式，沿用 /balance 的位置报错结构。
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

    # 骨架不得包含 H2O、H^+、e^-：它们由求解器作为带符号未知量补入。
    forbidden_errors = [
        {"side": side_name, "index": index, "formula": formula}
        for side_name, formulas in (
            ("reactants", request.reactants),
            ("products", request.products),
        )
        for index, formula in enumerate(formulas)
        if formula in ACIDIC_AUX_SPECIES
    ]
    if forbidden_errors:
        return JSONResponse(
            status_code=400,
            content={
                "balanced": False,
                "reason": "forbidden_species",
                "message": "半反应骨架不得包含 H2O、H^+、e^-，它们由求解器自动补入",
                "errors": forbidden_errors,
            },
        )

    try:
        result = balance_acidic_half_reaction(reactant, product)
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
    except ElectronCoefficientZeroError as exc:
        return JSONResponse(
            status_code=200,
            content={
                "balanced": False,
                "reason": "electron_coefficient_zero",
                "message": str(exc),
            },
        )

    return JSONResponse(status_code=200, content=result)


if __name__ == "__main__":
    import uvicorn

    host = os.environ.get("APP_HOST", "0.0.0.0")
    port = int(os.environ.get("APP_PORT", "8000"))
    uvicorn.run("app.main:app", host=host, port=port)
