/* Skaičiavimų modelis: koeficientai -> tikimybės, normaliojo skirstinio įvertis,
   rezultatų tinklelis ir sąlygų tikimybė. Be jokių puslapio (DOM) priklausomybių.
   Naršyklėje funkcijos tampa globaliomis, Node.js eksportuojamos per module.exports. */

// ---------- normalusis skirstinys ----------
function erf(x) { const s = Math.sign(x); x = Math.abs(x);
  const t = 1 / (1 + 0.3275911 * x);
  const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
  return s * y; }
const Phi = z => 0.5 * (1 + erf(z / Math.SQRT2));
function PhiInv(p) { // Acklam
  const a = [-39.69683028665376, 220.9460984245205, -275.9285104469687, 138.357751867269, -30.66479806614716, 2.506628277459239];
  const b = [-54.47609879822406, 161.5858368580409, -155.6989798598866, 66.80131188771972, -13.28068155288572];
  const c = [-0.007784894002430293, -0.3223964580411365, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783];
  const d = [0.007784695709041462, 0.3224671290700398, 2.445134137142996, 3.754408661907416];
  const pl = 0.02425;
  if (p < pl) { const q = Math.sqrt(-2 * Math.log(p)); return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1); }
  if (p > 1 - pl) { const q = Math.sqrt(-2 * Math.log(1 - p)); return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1); }
  const q = p - 0.5, r = q * q;
  return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
}

// Each row -> threshold k and P(X > k) without margin. X = home−away margin for hcp, points otherwise.
function rowCalc(m, r) {
  const l = parseFloat(String(r.l).replace(",", ".")), a = parseFloat(String(r.a).replace(",", ".")), b = parseFloat(String(r.b).replace(",", "."));
  if (!isFinite(l) || !(a > 1) || !(b > 1)) return null;
  const ia = 1 / a, ib = 1 / b, book = ia + ib;
  return { k: m.key === "hcp" ? -l : l, pa: ia / book, pb: ib / book, margin: book - 1 };
}
function fit(m, rows) {
  const pts = rows.map(r => rowCalc(m, r)).filter(Boolean).map(c => ({ k: c.k, z: PhiInv(1 - Math.min(0.99, Math.max(0.01, c.pa))) }));
  if (!pts.length) return null;
  if (pts.length >= 2) {
    const n = pts.length, mk = pts.reduce((s, p) => s + p.k, 0) / n, mz = pts.reduce((s, p) => s + p.z, 0) / n;
    const sxx = pts.reduce((s, p) => s + (p.k - mk) ** 2, 0);
    if (sxx > 0) {
      const slope = pts.reduce((s, p) => s + (p.k - mk) * (p.z - mz), 0) / sxx;
      if (slope > 0.01) { const sigma = 1 / slope; return { mu: mk - mz * sigma, sigma, n, assumed: false }; }
    }
  }
  const sigma = m.sigma0;
  const mu = pts.reduce((s, p) => s + (p.k - p.z * sigma), 0) / pts.length;
  return { mu, sigma, n: pts.length, assumed: true };
}
const pOver = (f, k) => 1 - Phi((k - f.mu) / f.sigma);

// ---------- rezultato modelis ----------
// F = {hcp, tot, home, away}: kiekvienos rinkos fit() rezultatas arba null.
function scoreModel(F) {
  if (F.hcp && F.tot) {
    const vT = F.tot.sigma ** 2, vM = F.hcp.sigma ** 2;
    const s = Math.sqrt((vT + vM) / 4);
    return { mh: (F.tot.mu + F.hcp.mu) / 2, ma: (F.tot.mu - F.hcp.mu) / 2, sh: s, sa: s, rho: (vT - vM) / (vT + vM), src: "handikapas + suma" };
  }
  if (F.home && F.away) return { mh: F.home.mu, ma: F.away.mu, sh: F.home.sigma, sa: F.away.sigma, rho: 0, src: "komandų taškai" };
  return null;
}

// Visi sveiki rezultatai h : a per ±5σ su tikimybėmis (suma = 1).
function gridCells(M) {
  const rho = Math.max(-0.95, Math.min(0.95, M.rho)), cells = []; let tot = 0;
  for (let h = Math.max(0, Math.floor(M.mh - 5 * M.sh)); h <= Math.ceil(M.mh + 5 * M.sh); h++)
    for (let a = Math.max(0, Math.floor(M.ma - 5 * M.sa)); a <= Math.ceil(M.ma + 5 * M.sa); a++) {
      const x = (h - M.mh) / M.sh, y = (a - M.ma) / M.sa;
      const d = Math.exp(-(x * x - 2 * rho * x * y + y * y) / (2 * (1 - rho * rho)));
      cells.push([h, a, d]); tot += d;
    }
  cells.forEach(c => c[2] /= tot);
  return cells;
}

// ---------- sąlygos ----------
const VARS = [["home", h => h], ["away", (h, a) => a], ["total", (h, a) => h + a], ["diff", (h, a) => h - a]];
const OPS = { ">": (x, y) => x > y, ">=": (x, y) => x >= y, "<": (x, y) => x < y, "<=": (x, y) => x <= y, "=": (x, y) => x === y };

// conds = [{v: "home"|"away"|"total"|"diff", op: ">"|">="|"<"|"<="|"=", x: skaičius}], visos sujungtos „ir“.
function conditionProbability(M, conds) {
  const cs = conds.map(c => ({ f: VARS.find(v => v[0] === c.v)?.[1], op: OPS[c.op], x: Number(c.x) }));
  if (cs.some(c => !c.f || !c.op || !isFinite(c.x))) throw new Error("Netinkama sąlyga");
  let p = 0;
  for (const [h, a, w] of gridCells(M)) if (cs.every(c => c.op(c.f(h, a), c.x))) p += w;
  return p;
}

if (typeof module !== "undefined") {
  module.exports = { erf, Phi, PhiInv, rowCalc, fit, pOver, scoreModel, gridCells, conditionProbability, VARS, OPS };
}
