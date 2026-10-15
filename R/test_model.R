# R modelio testai ir nepriklausomas Monte Carlo patikrinimas.
# Paleidimas iš repozitorijos šaknies:  Rscript R/test_model.R

source("R/model.R")

HCP <- matrix(c(-14.5, 2.45, 1.50,  -13.5, 2.30, 1.62,  -12.5, 2.15, 1.70,  -11.5, 1.95, 1.85,
                -10.5, 1.85, 1.95,   -9.5, 1.70, 2.15,   -8.5, 1.62, 2.30,   -7.5, 1.50, 2.45),
              ncol = 3, byrow = TRUE)
TOT <- matrix(c(160.5, 1.55, 2.35,  163.5, 1.75, 2.05,  165.5, 1.90, 1.90,
                167.5, 2.05, 1.75,  170.5, 2.35, 1.55), ncol = 3, byrow = TRUE)

near <- function(x, y, eps, msg) if (abs(x - y) >= eps) stop(sprintf("%s: %g != %g", msg, x, y))

r <- row_calc("hcp", -11.5, 1.95, 1.85)
near(r$pa, 0.4868, 1e-4, "P(namai -11.5)")
near(r$margin, 0.0534, 1e-4, "marža")

F <- list(hcp = fit_market("hcp", HCP), tot = fit_market("tot", TOT))
near(F$hcp$mu, 11.0, 0.05, "vid. skirtumas")
near(F$hcp$sigma, 11.31, 0.05, "σ skirtumo")
near(p_over(F$hcp, 0), 0.8346, 1e-3, "P(Fenerbahce laimi)")
near(F$tot$mu, 165.5, 0.05, "vid. suma")

M <- score_model(F)
near(M$mh, 88.25, 0.05, "namų vidurkis")
near(M$ma, 77.25, 0.05, "svečių vidurkis")
near(M$rho, 0.491, 0.01, "koreliacija")
near(sum(grid_cells(M)$p), 1, 1e-9, "tinklelio suma")

p_away <- condition_probability(M, list(list(v = "diff", op = "<", x = 0)))
p_home <- condition_probability(M, list(list(v = "diff", op = ">", x = 0)))
p_draw <- condition_probability(M, list(list(v = "diff", op = "=", x = 0)))
near(p_away + p_home + p_draw, 1, 1e-9, "laimi + pralaimi + lygiosios")
conds <- list(list(v = "home", op = ">", x = 75), list(v = "away", op = ">", x = 75), list(v = "diff", op = "<", x = 0))
p_both <- condition_probability(M, conds)
near(p_both, 0.1018, 2e-3, "abi > 75 ir laimi Partizan (sutampa su JS ir Python)")

# Nepriklausomas patikrinimas: 2 mln. dvimačio normaliojo imčių per Cholesky, suapvalinta iki sveikų taškų
set.seed(1)
n <- 2e6
z1 <- rnorm(n); z2 <- rnorm(n)
hs <- round(M$mh + M$sh * z1)
as <- round(M$ma + M$sa * (M$rho * z1 + sqrt(1 - M$rho^2) * z2))
mc <- mean(hs > 75 & as > 75 & hs < as)
near(p_both, mc, 1.5e-3, "sąlyga prieš Monte Carlo")
cat(sprintf("Monte Carlo palyginimas: sąlyga %.2f%% (modelis) vs %.2f%% (Monte Carlo)\n", p_both * 100, mc * 100))
cat("Visi R testai praėjo.\n")
