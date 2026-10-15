// Skaičiavimų patikra: node tests/model.test.js
const assert = require("node:assert/strict");
const M = require("../app/model.js");

const HCP = { key: "hcp", sigma0: 12 };
const TOT = { key: "tot", sigma0: 16 };
const rows = list => list.map(([l, a, b]) => ({ l, a, b }));

// Fenerbahce – Partizan, Betsson handikapai
const hcp = rows([[-14.5, 2.45, 1.50], [-13.5, 2.30, 1.62], [-12.5, 2.15, 1.70], [-11.5, 1.95, 1.85],
                  [-10.5, 1.85, 1.95], [-9.5, 1.70, 2.15], [-8.5, 1.62, 2.30], [-7.5, 1.50, 2.45]]);
const tot = rows([[160.5, 1.55, 2.35], [163.5, 1.75, 2.05], [165.5, 1.90, 1.90], [167.5, 2.05, 1.75], [170.5, 2.35, 1.55]]);

const near = (x, y, eps, msg) => assert.ok(Math.abs(x - y) < eps, `${msg}: ${x} != ${y}`);

// Normalusis skirstinys
near(M.Phi(1.96), 0.975, 1e-4, "Phi(1.96)");
near(M.PhiInv(0.975), 1.96, 1e-3, "PhiInv(0.975)");

// Koeficientai -> tikimybės be maržos
const r = M.rowCalc(HCP, hcp[3]);
near(r.pa, 0.4868, 1e-4, "P(namai -11.5)");
near(r.margin, 0.0534, 1e-4, "marža");

// Normaliojo skirstinio pritaikymas
const F = { hcp: M.fit(HCP, hcp), tot: M.fit(TOT, tot), home: null, away: null };
near(F.hcp.mu, 11.0, 0.05, "vid. skirtumas");
near(F.hcp.sigma, 11.31, 0.05, "σ skirtumo");
near(M.pOver(F.hcp, 0), 0.8346, 1e-3, "P(Fenerbahce laimi)");
near(F.tot.mu, 165.5, 0.05, "vid. suma");

// Rezultato modelis ir tinklelis
const SM = M.scoreModel(F);
near(SM.mh, 88.25, 0.05, "namų vidurkis");
near(SM.ma, 77.25, 0.05, "svečių vidurkis");
near(SM.rho, 0.491, 0.01, "koreliacija");
const cells = M.gridCells(SM);
near(cells.reduce((s, c) => s + c[2], 0), 1, 1e-9, "tinklelio suma");

// Sąlygos
const pAway = M.conditionProbability(SM, [{ v: "diff", op: "<", x: 0 }]);
const pHome = M.conditionProbability(SM, [{ v: "diff", op: ">", x: 0 }]);
const pDraw = M.conditionProbability(SM, [{ v: "diff", op: "=", x: 0 }]);
near(pAway + pHome + pDraw, 1, 1e-9, "laimi + pralaimi + lygiosios");
const pBoth = M.conditionProbability(SM, [{ v: "home", op: ">", x: 75 }, { v: "away", op: ">", x: 75 }, { v: "diff", op: "<", x: 0 }]);
assert.ok(pBoth < pAway, "papildomos sąlygos mažina tikimybę");
near(pBoth, 0.1018, 2e-3, "abi > 75 ir laimi Partizan");
assert.throws(() => M.conditionProbability(SM, [{ v: "nėra", op: ">", x: 1 }]));

// CSV įkėlimas
const fs = require("node:fs");
const path = require("node:path");
const own = M.parseOddsCSV(fs.readFileSync(path.join(__dirname, "../data/fenerbahce_partizan_betsson_handicap.csv"), "utf8"));
assert.equal(own.count, 8, "esamas handikapo failas");
assert.deepEqual(own.byMarket.hcp?.[3], { l: "-11.5", a: "1.95", b: "1.85" }, "market=Handicap -> hcp");

const semi = M.parseOddsCSV("Linija;Daugiau;Mažiau\n165,5;1,90;1,90\n167,5;2,05;1,75\n");
assert.equal(semi.unassigned.length, 2, "kabliataškis ir dešimtainis kablelis");
assert.deepEqual(semi.unassigned[0], { l: "165.5", a: "1.90", b: "1.90" });

const combo = M.parseOddsCSV(["rinka,linija,a,b", "handikapas,-11.5,1.95,1.85", "taškų suma,165.5,1.90,1.90",
  "namai,88.5,1.85,1.95", "svečiai,77.5,1.90,1.90", "tekstas be skaičių"].join("\n"));
assert.equal(combo.byMarket.hcp.length, 1); assert.equal(combo.byMarket.tot.length, 1);
assert.equal(combo.byMarket.home.length, 1); assert.equal(combo.byMarket.away.length, 1);
assert.equal(combo.count, 4, "eilutės be skaičių praleidžiamos");

console.log("Visi skaičiavimų testai praėjo.");
