"""Skaičiavimų modelis Python kalba (tas pats kaip app/model.js).

Koeficientai -> tikimybės be maržos, normaliojo skirstinio įvertis kiekvienai rinkai,
abiejų komandų taškų modelis su koreliacija, rezultatų tinklelis ir sąlygų tikimybė.
Naudoja tik Python standartinę biblioteką.

Paleidimas iš repozitorijos šaknies:
    python3 python/model.py --hcp data/hcp.csv --tot data/tot.csv \
        --home-name Fenerbahce --away-name Partizan \
        --cond "home>75" --cond "away>75" --cond "diff<0"

CSV failuose naudojami pirmi trys skaitiniai stulpeliai: linija, koeficientas A, koeficientas B.
Handikapui linija yra namų komandos handikapas (pvz. -11.5), A = namai, B = svečiai.
Sumoms ir komandų taškams A = daugiau, B = mažiau.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from statistics import NormalDist

N = NormalDist()
SIGMA0 = {"hcp": 12.0, "tot": 16.0, "home": 10.0, "away": 10.0}


# ---------- koeficientai -> tikimybės ----------
@dataclass
class Row:
    k: float        # slenkstis: X > k (handikapui X = namai − svečiai, k = −linija)
    pa: float       # P(A) be maržos
    pb: float       # P(B) be maržos
    margin: float   # bukmekerio marža


def row_calc(market: str, line: float, a: float, b: float) -> Row | None:
    if not (a > 1 and b > 1) or not math.isfinite(line):
        return None
    ia, ib = 1 / a, 1 / b
    book = ia + ib
    return Row(k=-line if market == "hcp" else line, pa=ia / book, pb=ib / book, margin=book - 1)


# ---------- normaliojo skirstinio įvertis ----------
@dataclass
class Fit:
    mu: float
    sigma: float
    n: int
    assumed: bool   # True, kai σ neįvertinta (viena linija) ir paimta numatyta


def fit(market: str, rows: list[tuple[float, float, float]]) -> Fit | None:
    """P(X > k) = 1 − Φ((k − μ)/σ): z = Φ⁻¹(1 − p) tiesiškai priklauso nuo k, nuolydis = 1/σ."""
    calc = [c for c in (row_calc(market, *r) for r in rows) if c]
    if not calc:
        return None
    pts = [(c.k, N.inv_cdf(1 - min(0.99, max(0.01, c.pa)))) for c in calc]
    if len(pts) >= 2:
        n = len(pts)
        mk = sum(k for k, _ in pts) / n
        mz = sum(z for _, z in pts) / n
        sxx = sum((k - mk) ** 2 for k, _ in pts)
        if sxx > 0:
            slope = sum((k - mk) * (z - mz) for k, z in pts) / sxx
            if slope > 0.01:
                sigma = 1 / slope
                return Fit(mu=mk - mz * sigma, sigma=sigma, n=n, assumed=False)
    sigma = SIGMA0[market]
    mu = sum(k - z * sigma for k, z in pts) / len(pts)
    return Fit(mu=mu, sigma=sigma, n=len(pts), assumed=True)


def p_over(f: Fit, k: float) -> float:
    return 1 - N.cdf((k - f.mu) / f.sigma)


# ---------- rezultato modelis ----------
@dataclass
class ScoreModel:
    mh: float   # namų taškų vidurkis
    ma: float   # svečių taškų vidurkis
    sh: float   # namų σ
    sa: float   # svečių σ
    rho: float  # koreliacija tarp komandų taškų
    src: str


def score_model(F: dict[str, Fit | None]) -> ScoreModel | None:
    """Iš handikapo ir sumos: var(H) = var(A) = (σT² + σM²)/4, ρ = (σT² − σM²)/(σT² + σM²)."""
    hcp, tot, home, away = (F.get(k) for k in ("hcp", "tot", "home", "away"))
    if hcp and tot:
        vT, vM = tot.sigma ** 2, hcp.sigma ** 2
        s = math.sqrt((vT + vM) / 4)
        return ScoreModel(mh=(tot.mu + hcp.mu) / 2, ma=(tot.mu - hcp.mu) / 2, sh=s, sa=s,
                          rho=(vT - vM) / (vT + vM), src="handikapas + suma")
    if home and away:
        return ScoreModel(mh=home.mu, ma=away.mu, sh=home.sigma, sa=away.sigma, rho=0.0, src="komandų taškai")
    return None


def grid_cells(M: ScoreModel) -> list[tuple[int, int, float]]:
    """Visi sveiki rezultatai h : a per ±5σ su dvimačio normaliojo tankiu, sunormuoti iki 1."""
    rho = max(-0.95, min(0.95, M.rho))
    cells, tot = [], 0.0
    for h in range(max(0, math.floor(M.mh - 5 * M.sh)), math.ceil(M.mh + 5 * M.sh) + 1):
        for a in range(max(0, math.floor(M.ma - 5 * M.sa)), math.ceil(M.ma + 5 * M.sa) + 1):
            x, y = (h - M.mh) / M.sh, (a - M.ma) / M.sa
            d = math.exp(-(x * x - 2 * rho * x * y + y * y) / (2 * (1 - rho * rho)))
            cells.append((h, a, d))
            tot += d
    return [(h, a, d / tot) for h, a, d in cells]


# ---------- sąlygos ----------
VARS = {"home": lambda h, a: h, "away": lambda h, a: a, "total": lambda h, a: h + a, "diff": lambda h, a: h - a}
OPS = {">": lambda x, y: x > y, ">=": lambda x, y: x >= y, "<": lambda x, y: x < y,
       "<=": lambda x, y: x <= y, "=": lambda x, y: x == y}


def condition_probability(M: ScoreModel, conds: list[tuple[str, str, float]]) -> float:
    """conds = [("home", ">", 75), ...], visos sujungtos „ir“."""
    for v, op, _ in conds:
        if v not in VARS or op not in OPS:
            raise ValueError(f"Netinkama sąlyga: {v} {op}")
    return sum(w for h, a, w in grid_cells(M)
               if all(OPS[op](VARS[v](h, a), x) for v, op, x in conds))


def parse_cond(text: str) -> tuple[str, str, float]:
    m = re.fullmatch(r"\s*(home|away|total|diff)\s*(>=|<=|>|<|=)\s*(-?\d+(?:[.,]\d+)?)\s*", text)
    if not m:
        raise argparse.ArgumentTypeError(f"Sąlyga turi būti pvz. home>75, diff<0: {text!r}")
    return m.group(1), m.group(2), float(m.group(3).replace(",", "."))


# ---------- CSV ir komandinė eilutė ----------
def read_rows(path: str) -> list[tuple[float, float, float]]:
    out = []
    with open(path, encoding="utf-8-sig", newline="") as f:
        for rec in csv.reader(f):
            nums = []
            for cell in rec:
                try:
                    nums.append(float(cell.replace(",", ".")))
                except ValueError:
                    pass
            if len(nums) >= 3:
                out.append((nums[0], nums[1], nums[2]))
    return out


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description="Rungtynių tikimybės iš lažybų koeficientų")
    for k, name in (("hcp", "handikapas"), ("tot", "taškų suma"), ("home", "namų taškai"), ("away", "svečių taškai")):
        ap.add_argument(f"--{k}", metavar="CSV", help=f"{name}: linija, A, B")
    ap.add_argument("--home-name", default="Namai")
    ap.add_argument("--away-name", default="Svečiai")
    ap.add_argument("--cond", type=parse_cond, action="append", default=[],
                    help="sąlyga, pvz. home>75 (galima kartoti, visos sujungtos „ir“)")
    ap.add_argument("--top", type=int, default=10, help="kiek tikėtiniausių rezultatų rodyti")
    args = ap.parse_args(argv)

    F = {k: fit(k, read_rows(getattr(args, k))) if getattr(args, k) else None for k in ("hcp", "tot", "home", "away")}
    H, A = args.home_name, args.away_name
    if F["hcp"]:
        print(f"{H} laimi: {p_over(F['hcp'], 0) * 100:.2f}%   skirtumas {F['hcp'].mu:+.1f} (σ {F['hcp'].sigma:.1f})")
    if F["tot"]:
        print(f"Taškų suma: {F['tot'].mu:.1f} (σ {F['tot'].sigma:.1f})")
    M = score_model(F)
    if not M:
        print("Rezultato modeliui reikia handikapo ir sumos (arba abiejų komandų taškų).")
        return
    print(f"Tikėtinas rezultatas: {H} {M.mh:.1f} : {M.ma:.1f} {A}   ρ = {M.rho:.3f}   ({M.src})")
    cells = sorted((c for c in grid_cells(M) if c[0] != c[1]), key=lambda c: -c[2])[: args.top]
    print("\nTikėtiniausi rezultatai:")
    for i, (h, a, p) in enumerate(cells, 1):
        print(f"  {i:>2}. {h} : {a}   {p * 100:.2f}%")
    if args.cond:
        p = condition_probability(M, args.cond)
        conds = " ir ".join(f"{v} {op} {x:g}" for v, op, x in args.cond)
        print(f"\nP({conds}) = {p * 100:.2f}%   teisingas koeficientas {1 / p:.2f}" if p > 0 else f"\nP({conds}) = 0")


if __name__ == "__main__":
    main()
