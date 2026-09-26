"""POST /balance/acidic-half 的接口测试。"""

from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)

PATH = "/balance/acidic-half"


def test_acidic_half_success_reduction():
    response = client.post(
        PATH, json={"reactants": ["MnO4^-"], "products": ["Mn^2+"]}
    )
    assert response.status_code == 200
    data = response.json()
    assert data["balanced"] is True
    assert data["coefficients"] == {"reactants": [1, 8, 5], "products": [1, 4]}
    assert data["equation"] == "MnO4^- + 8H^+ + 5e^- -> Mn^2+ + 4H2O"
    assert data["atom_totals"]["reactants"] == {"H": 8, "Mn": 1, "O": 4}
    assert data["atom_totals"]["products"] == {"H": 8, "Mn": 1, "O": 4}
    assert data["charge_totals"] == {"reactants": 2, "products": 2}


def test_acidic_half_success_oxidation_auxiliaries_sorted():
    # 氧化半反应：水落在反应物侧，H^+ 与 e^- 落在生成物侧
    response = client.post(
        PATH, json={"reactants": ["NH4^+"], "products": ["NO3^-"]}
    )
    assert response.status_code == 200
    data = response.json()
    assert data["balanced"] is True
    assert [s["formula"] for s in data["reactants"]] == ["NH4^+", "H2O"]
    assert [s["formula"] for s in data["products"]] == ["NO3^-", "H^+", "e^-"]
    assert data["equation"] == "NH4^+ + 3H2O -> NO3^- + 10H^+ + 8e^-"


def test_acidic_half_simple_electron_only():
    response = client.post(
        PATH, json={"reactants": ["Fe^2+"], "products": ["Fe^3+"]}
    )
    assert response.status_code == 200
    data = response.json()
    assert data["balanced"] is True
    assert data["coefficients"] == {"reactants": [1], "products": [1, 1]}
    assert data["equation"] == "Fe^2+ -> Fe^3+ + e^-"
    assert data["charge_totals"] == {"reactants": 2, "products": 2}


def test_acidic_half_electron_coefficient_zero():
    response = client.post(
        PATH, json={"reactants": ["P4O10"], "products": ["H3PO4"]}
    )
    assert response.status_code == 200
    data = response.json()
    assert data["balanced"] is False
    assert data["reason"] == "electron_coefficient_zero"


def test_acidic_half_multiple_solutions():
    response = client.post(
        PATH, json={"reactants": ["H2O2"], "products": ["O2"]}
    )
    assert response.status_code == 200
    data = response.json()
    assert data["balanced"] is False
    assert data["reason"] == "multiple_solutions"


def test_acidic_half_unbalanceable():
    response = client.post(
        PATH, json={"reactants": ["NaCl"], "products": ["Cl^-"]}
    )
    assert response.status_code == 200
    data = response.json()
    assert data["balanced"] is False
    assert data["reason"] == "unbalanceable"


def test_acidic_half_forbidden_water_in_skeleton():
    response = client.post(PATH, json={"reactants": ["H2O"], "products": ["H2"]})
    assert response.status_code == 400
    data = response.json()
    assert data["balanced"] is False
    assert data["reason"] == "forbidden_species"
    assert data["errors"] == [
        {"side": "reactants", "index": 0, "formula": "H2O"}
    ]


def test_acidic_half_forbidden_h_plus_and_electron():
    response = client.post(
        PATH, json={"reactants": ["H^+"], "products": ["e^-"]}
    )
    assert response.status_code == 400
    data = response.json()
    assert data["reason"] == "forbidden_species"
    assert {(e["side"], e["formula"]) for e in data["errors"]} == {
        ("reactants", "H^+"),
        ("products", "e^-"),
    }


def test_acidic_half_invalid_formula_reports_position():
    response = client.post(PATH, json={"reactants": ["Na+"], "products": ["Na"]})
    assert response.status_code == 400
    data = response.json()
    assert data["reason"] == "invalid_formula"
    assert data["errors"][0]["position"] == 2


def test_acidic_half_formula_check_takes_precedence_over_forbidden():
    # 同一请求里不会出现既非法又是禁词的化学式；非法化学式先报 invalid_formula
    response = client.post(PATH, json={"reactants": ["X2"], "products": ["e^-"]})
    assert response.status_code == 400
    data = response.json()
    assert data["reason"] == "invalid_formula"


def test_acidic_half_two_species_per_side_returns_422():
    response = client.post(
        PATH, json={"reactants": ["Fe", "H2"], "products": ["Fe^3+"]}
    )
    assert response.status_code == 422


def test_acidic_half_empty_side_returns_422():
    response = client.post(PATH, json={"reactants": [], "products": ["Fe^2+"]})
    assert response.status_code == 422


def test_acidic_half_missing_field_returns_422():
    response = client.post(PATH, json={"reactants": ["Fe"]})
    assert response.status_code == 422


def test_original_balance_endpoint_unchanged():
    # 回归：原有 /balance 行为不受影响
    response = client.post(
        "/balance",
        json={"reactants": ["Fe", "O2"], "products": ["Fe3O4"]},
    )
    assert response.status_code == 200
    assert response.json()["coefficients"] == {
        "reactants": [3, 2],
        "products": [1],
    }
