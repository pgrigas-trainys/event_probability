"""Python modelio testai ir palyginimas su SciPy.

Paleidimas iš repozitorijos šaknies:  python3 python/test_model.py
SciPy dalis praleidžiama, jei SciPy neįdiegta (pip install scipy).
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import model as m  # noqa: E402

HCP = [(-14.5, 2.45, 1.50), (-13.5, 2.30, 1.62), (-12.5, 2.15, 1.70), (-11.5, 1.95, 1.85),
       (-10.5, 1.85, 1.95), (-9.5, 1.70, 2.15), (-8.5, 1.62, 2.30), (-7.5, 1.50, 2.45)]
TOT = [(160.5, 1.55, 2.35), (163.5, 1.75, 2.05), (165.5, 1.90, 1.90), (167.5, 2.05, 1.75), (170.5, 2.35, 1.55)]


def near(x, y, eps, msg):
    assert abs(x - y) < eps, f"{msg}: {x} != {y}"


def test_model():
    r = m.row_calc("hcp", -11.5, 1.95, 1.85)
    near(r.pa, 0.4868, 1e-4, "P(namai -11.5)")
    near(r.margin, 0.0534, 1e-4, "marža")

    F = {"hcp": m.fit("hcp", HCP), "tot": m.fit("tot", TOT), "home": None, "away": None}
    near(F["hcp"].mu, 11.0, 0.05, "vid. skirtumas")
    near(F["hcp"].sigma, 11.31, 0.05, "σ skirtumo")
    near(m.p_over(F["hcp"], 0), 0.8346, 1e-3, "P(Fenerbahce laimi)")
    near(F["tot"].mu, 165.5, 0.05, "vid. suma")

    M = m.score_model(F)
    near(M.mh, 88.25, 0.05, "namų vidurkis")
    near(M.ma, 77.25, 0.05, "svečių vidurkis")
    near(M.rho, 0.491, 0.01, "koreliacija")
    near(sum(c[2] for c in m.grid_cells(M)), 1, 1e-9, "tinklelio suma")

    p_away = m.condition_probability(M, [("diff", "<", 0)])
    p_home = m.condition_probability(M, [("diff", ">", 0)])
    p_draw = m.condition_probability(M, [("diff", "=", 0)])
    near(p_away + p_home + p_draw, 1, 1e-9, "laimi + pralaimi + lygiosios")
    p_both = m.condition_probability(M, [("home", ">", 75), ("away", ">", 75), ("diff", "<", 0)])
    near(p_both, 0.1018, 2e-3, "abi > 75 ir laimi Partizan (sutampa su JS)")
    return F, M


def test_against_scipy(F, M):
    try:
        import numpy as np
        from scipy import stats
    except ImportError:
        print("SciPy neįdiegta: palyginimas praleistas.")
        return
    # Φ ir Φ⁻¹
    near(m.N.cdf(1.2345), stats.norm.cdf(1.2345), 1e-9, "Φ")
    near(m.N.inv_cdf(0.3), stats.norm.ppf(0.3), 1e-9, "Φ⁻¹")

    # Rezultato tinklelis prieš SciPy dvimatį normalųjį: tikimybė langeliui h : a
    # integruojama per [h−0.5, h+0.5] × [a−0.5, a+0.5].
    cov = [[M.sh ** 2, M.rho * M.sh * M.sa], [M.rho * M.sh * M.sa, M.sa ** 2]]
    mvn = stats.multivariate_normal(mean=[M.mh, M.ma], cov=cov)

    def box(h0, h1, a0, a1):
        return mvn.cdf([h1, a1]) - mvn.cdf([h0, a1]) - mvn.cdf([h1, a0]) + mvn.cdf([h0, a0])

    grid = {(h, a): p for h, a, p in m.grid_cells(M)}
    for h, a in [(88, 77), (95, 70), (80, 82)]:
        near(grid[(h, a)], box(h - 0.5, h + 0.5, a - 0.5, a + 0.5), 2e-5, f"langelis {h}:{a}")

    # Sąlyga „abi > 75 ir laimi svečiai“ Monte Carlo su SciPy (2 mln. imčių, suapvalinta iki sveikų taškų)
    rng = np.random.default_rng(1)
    xs = np.rint(mvn.rvs(size=2_000_000, random_state=rng))
    mc = np.mean((xs[:, 0] > 75) & (xs[:, 1] > 75) & (xs[:, 0] < xs[:, 1]))
    p = m.condition_probability(M, [("home", ">", 75), ("away", ">", 75), ("diff", "<", 0)])
    near(p, mc, 1.5e-3, "sąlyga prieš Monte Carlo")
    print(f"SciPy palyginimas: sąlyga {p * 100:.2f}% (modelis) vs {mc * 100:.2f}% (Monte Carlo)")


if __name__ == "__main__":
    F, M = test_model()
    test_against_scipy(F, M)
    print("Visi Python testai praėjo.")
