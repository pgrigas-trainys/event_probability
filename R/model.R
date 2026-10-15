# Skaičiavimų modelis R kalba (tas pats kaip app/model.js ir python/model.py).
#
# Koeficientai -> tikimybės be maržos, normaliojo skirstinio įvertis kiekvienai rinkai,
# abiejų komandų taškų modelis su koreliacija, rezultatų tinklelis ir sąlygų tikimybė.
# Naudojamas tik bazinis R.
#
# Paleidimas iš repozitorijos šaknies:
#   Rscript R/model.R --hcp data/hcp.csv --tot data/tot.csv \
#     --home-name Fenerbahce --away-name Partizan \
#     --cond "home>75" --cond "away>75" --cond "diff<0"
#
# Naudojimas iš kito R kodo:  source("R/model.R")
#
# CSV failuose naudojami pirmi trys skaitiniai stulpeliai: linija, koeficientas A, koeficientas B.
# Handikapui linija yra namų komandos handikapas (pvz. -11.5), A = namai, B = svečiai.
# Sumoms ir komandų taškams A = daugiau, B = mažiau.

SIGMA0 <- c(hcp = 12, tot = 16, home = 10, away = 10)

# ---------- koeficientai -> tikimybės ----------
row_calc <- function(market, line, a, b) {
  if (!is.finite(line) || !(a > 1 && b > 1)) return(NULL)
  ia <- 1 / a; ib <- 1 / b; book <- ia + ib
  list(k = if (market == "hcp") -line else line,   # slenkstis: X > k
       pa = ia / book, pb = ib / book, margin = book - 1)
}

# ---------- normaliojo skirstinio įvertis ----------
# P(X > k) = 1 - Φ((k - μ)/σ): z = Φ⁻¹(1 - p) tiesiškai priklauso nuo k, nuolydis = 1/σ.
fit_market <- function(market, rows) {
  calc <- Filter(Negate(is.null), lapply(seq_len(nrow(rows)), function(i)
    row_calc(market, rows[i, 1], rows[i, 2], rows[i, 3])))
  if (!length(calc)) return(NULL)
  k <- vapply(calc, `[[`, 0, "k")
  z <- qnorm(1 - pmin(0.99, pmax(0.01, vapply(calc, `[[`, 0, "pa"))))
  n <- length(k)
  if (n >= 2) {
    sxx <- sum((k - mean(k))^2)
    if (sxx > 0) {
      slope <- sum((k - mean(k)) * (z - mean(z))) / sxx
      if (slope > 0.01) {
        sigma <- 1 / slope
        return(list(mu = mean(k) - mean(z) * sigma, sigma = sigma, n = n, assumed = FALSE))
      }
    }
  }
  sigma <- SIGMA0[[market]]
  list(mu = mean(k - z * sigma), sigma = sigma, n = n, assumed = TRUE)
}

p_over <- function(f, k) 1 - pnorm((k - f$mu) / f$sigma)

# ---------- rezultato modelis ----------
# Iš handikapo ir sumos: var(H) = var(A) = (σT² + σM²)/4, ρ = (σT² − σM²)/(σT² + σM²).
score_model <- function(F) {
  if (!is.null(F$hcp) && !is.null(F$tot)) {
    vT <- F$tot$sigma^2; vM <- F$hcp$sigma^2; s <- sqrt((vT + vM) / 4)
    return(list(mh = (F$tot$mu + F$hcp$mu) / 2, ma = (F$tot$mu - F$hcp$mu) / 2,
                sh = s, sa = s, rho = (vT - vM) / (vT + vM), src = "handikapas + suma"))
  }
  if (!is.null(F$home) && !is.null(F$away))
    return(list(mh = F$home$mu, ma = F$away$mu, sh = F$home$sigma, sa = F$away$sigma,
                rho = 0, src = "komandų taškai"))
  NULL
}

# Visi sveiki rezultatai h : a per ±5σ su dvimačio normaliojo tankiu, sunormuoti iki 1.
grid_cells <- function(M) {
  rho <- max(-0.95, min(0.95, M$rho))
  h <- seq(max(0, floor(M$mh - 5 * M$sh)), ceiling(M$mh + 5 * M$sh))
  a <- seq(max(0, floor(M$ma - 5 * M$sa)), ceiling(M$ma + 5 * M$sa))
  g <- expand.grid(h = h, a = a)
  x <- (g$h - M$mh) / M$sh; y <- (g$a - M$ma) / M$sa
  d <- exp(-(x * x - 2 * rho * x * y + y * y) / (2 * (1 - rho * rho)))
  g$p <- d / sum(d)
  g
}

# ---------- sąlygos ----------
VARS <- list(home = function(h, a) h, away = function(h, a) a,
             total = function(h, a) h + a, diff = function(h, a) h - a)
OPS <- list(">" = `>`, ">=" = `>=`, "<" = `<`, "<=" = `<=`, "=" = `==`)

# conds: sąrašas list(v = "home", op = ">", x = 75), visos sujungtos „ir“.
condition_probability <- function(M, conds) {
  g <- grid_cells(M)
  ok <- rep(TRUE, nrow(g))
  for (c in conds) {
    if (is.null(VARS[[c$v]]) || is.null(OPS[[c$op]])) stop("Netinkama sąlyga: ", c$v, " ", c$op)
    ok <- ok & OPS[[c$op]](VARS[[c$v]](g$h, g$a), c$x)
  }
  sum(g$p[ok])
}

parse_cond <- function(text) {
  m <- regmatches(text, regexec("^\\s*(home|away|total|diff)\\s*(>=|<=|>|<|=)\\s*(-?[0-9]+(?:[.,][0-9]+)?)\\s*$", text))[[1]]
  if (!length(m)) stop("Sąlyga turi būti pvz. home>75, diff<0: ", text)
  list(v = m[2], op = m[3], x = as.numeric(sub(",", ".", m[4])))
}

# ---------- CSV ir komandinė eilutė ----------
read_rows <- function(path) {
  lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
  out <- lapply(strsplit(lines, ","), function(cells) {
    nums <- suppressWarnings(as.numeric(sub(",", ".", trimws(cells))))
    nums <- nums[!is.na(nums)]
    if (length(nums) >= 3) nums[1:3] else NULL
  })
  out <- Filter(Negate(is.null), out)
  if (!length(out)) return(matrix(numeric(0), ncol = 3))
  do.call(rbind, out)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  opt <- list(cond = list(), home_name = "Namai", away_name = "Svečiai", top = 10)
  i <- 1
  while (i <= length(args)) {
    key <- args[i]; val <- args[i + 1]
    switch(key,
      "--hcp" = , "--tot" = , "--home" = , "--away" = { opt[[sub("--", "", key)]] <- val },
      "--home-name" = { opt$home_name <- val },
      "--away-name" = { opt$away_name <- val },
      "--cond" = { opt$cond[[length(opt$cond) + 1]] <- parse_cond(val) },
      "--top" = { opt$top <- as.integer(val) },
      stop("Nežinomas parametras: ", key))
    i <- i + 2
  }
  F <- list()
  for (k in c("hcp", "tot", "home", "away"))
    if (!is.null(opt[[k]])) F[k] <- list(fit_market(k, read_rows(opt[[k]])))
  H <- opt$home_name; A <- opt$away_name
  if (!is.null(F$hcp))
    cat(sprintf("%s laimi: %.2f%%   skirtumas %+.1f (σ %.1f)\n", H, p_over(F$hcp, 0) * 100, F$hcp$mu, F$hcp$sigma))
  if (!is.null(F$tot)) cat(sprintf("Taškų suma: %.1f (σ %.1f)\n", F$tot$mu, F$tot$sigma))
  M <- score_model(F)
  if (is.null(M)) { cat("Rezultato modeliui reikia handikapo ir sumos (arba abiejų komandų taškų).\n"); return(invisible()) }
  cat(sprintf("Tikėtinas rezultatas: %s %.1f : %.1f %s   ρ = %.3f   (%s)\n", H, M$mh, M$ma, A, M$rho, M$src))
  g <- grid_cells(M); g <- g[g$h != g$a, ]; g <- head(g[order(-g$p), ], opt$top)
  cat("\nTikėtiniausi rezultatai:\n")
  for (j in seq_len(nrow(g))) cat(sprintf("  %2d. %d : %d   %.2f%%\n", j, g$h[j], g$a[j], g$p[j] * 100))
  if (length(opt$cond)) {
    p <- condition_probability(M, opt$cond)
    txt <- paste(vapply(opt$cond, function(c) sprintf("%s %s %g", c$v, c$op, c$x), ""), collapse = " ir ")
    if (p > 0) cat(sprintf("\nP(%s) = %.2f%%   teisingas koeficientas %.2f\n", txt, p * 100, 1 / p))
    else cat(sprintf("\nP(%s) = 0\n", txt))
  }
}

if (sys.nframe() == 0) main()
