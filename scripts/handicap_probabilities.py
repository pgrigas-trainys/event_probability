"""Handikapų koeficientai -> tikimybės be maržos + normaliojo skirstinio įvertis.

Paleidimas: python3 scripts/handicap_probabilities.py data/<failas>.csv
"""
import csv
import math
import sys
from statistics import NormalDist


def main(path):
    with open(path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    out = []
    for r in rows:
        line = float(r["line_fenerbahce"])
        o1, o2 = float(r["odds_fenerbahce"]), float(r["odds_partizan"])
        i1, i2 = 1 / o1, 1 / o2
        book = i1 + i2
        out.append({"k": -line, "o1": o1, "o2": o2, "p": i1 / book, "margin": book - 1})

    # Normal fit: P(margin > k) = 1 - Phi((k - mu) / sigma); linear in z.
    xs = [r["k"] for r in out]
    zs = [NormalDist().inv_cdf(1 - r["p"]) for r in out]  # z = (k - mu)/sigma
    n = len(xs); mx = sum(xs) / n; mz = sum(zs) / n
    b = sum((x - mx) * (z - mz) for x, z in zip(xs, zs)) / sum((x - mx) ** 2 for x in xs)
    a = mz - b * mx
    sigma = 1 / b; mu = -a / b
    nd = NormalDist(mu, sigma)

    print(f"{'Fener laimi':<22}{'Koef.':>7}{'Tikimybė':>10}{'Marža':>8}")
    for r in out:
        print(f"{'daugiau nei ' + str(r['k']) + ' tšk.':<22}{r['o1']:>7.2f}{r['p']*100:>9.1f}%{r['margin']*100:>7.1f}%")
    print(f"\nĮvertis (normalusis skirstinys): Fenerbahce pranašumas vid. {mu:.1f} tšk., σ = {sigma:.1f}")
    print(f"P(Fenerbahce laimi) ≈ {(1 - nd.cdf(0)) * 100:.1f}%")
    for k in (0, 5, 10, 15, 20):
        print(f"  P(laimi > {k} tšk.) ≈ {(1 - nd.cdf(k)) * 100:.1f}%")


if __name__ == "__main__":
    main(sys.argv[1])
