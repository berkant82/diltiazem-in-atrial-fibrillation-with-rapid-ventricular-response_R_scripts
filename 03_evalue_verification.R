## =====================================================================
## DILTIAZEM FOR ATRIAL FIBRILLATION WITH RAPID VENTRICULAR RESPONSE
## Script 3 of 3  --  E-VALUES, WITH AN INDEPENDENT RE-DERIVATION
## =====================================================================
## This script deliberately does NOT reuse anything built by script 01.
## It reads the raw data file again, rebuilds the nine covariates and the
## outcome from the raw columns, and checks every derived quantity against
## the values printed in analysis_log_v2.txt before it reports anything.
## That is the point of the script: it is an independent reproduction of
## the weighted risk ratios and their E-values, not a shortcut.
##
## It produces ONLY the stabilised-IPTW risk ratios (log-link Poisson,
## robust variance) and the corresponding E-values for the three pairwise
## comparisons. Nothing else is recomputed.
##
## Read the VERDICT line at the end before using any number from here.
##
## HOW TO RUN
##   Put this script in the same folder as (or in a sub-folder next to)
##   the data file "diltiazem_doz1set_eng_.csv" and run it. Nothing has
##   to be edited.
##
##   From a terminal:   Rscript 03_evalue_verification.R
##   From RStudio:      open the file and press "Source"
##
## OUTPUT
##   output_v2/e_value.log   (plain ASCII)
##
## Written for R >= 4.1. Tested with R 4.3.
## =====================================================================


## ---------------------------------------------------------------------
## 0. SETTINGS -- normally nothing here needs to be edited
## ---------------------------------------------------------------------
DATA_NAME <- "diltiazem_doz1set_eng_.csv"
LOGFILE   <- "e_value.log"

## Expected values, copied from analysis_log_v2.txt. Every one is checked.
EXP <- list(
  n        = c(105, 59, 50),
  success  = c(58, 53, 44),
  male     = c(35, 26, 23),
  cad      = c(9, 5, 7),
  newaf    = c(27, 12, 14),
  ratectl  = c(80, 48, 37),
  age      = c(65.65, 62.22, 66.38),
  weight   = c(74.57, 75.53, 75.42),
  pulse    = c(149.22, 148.61, 153.04),
  sbp      = c(130.19, 131.36, 131.20),
  map      = c(98.57, 99.49, 98.93)
)
REF <- list(ld_sb = c(est = 0.6168009, lo = 0.5071233, hi = 0.7501989),
            ld_si = c(est = 0.6303746, lo = 0.5126439, hi = 0.7751427))


## ---------------------------------------------------------------------
## 1. ENVIRONMENT
## ---------------------------------------------------------------------
options(stringsAsFactors = FALSE, warn = 1, OutDec = ".")
invisible(suppressWarnings(Sys.setlocale("LC_COLLATE", "C")))   # locale-independent sorting
Sys.setenv(LANGUAGE = "en")

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit  <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(hit)) return(normalizePath(dirname(hit[1]), winslash = "/", mustWork = FALSE))
  for (i in rev(seq_len(sys.nframe()))) {
    of <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (is.character(of) && length(of) == 1L && nzchar(of))
      return(normalizePath(dirname(of), winslash = "/", mustWork = FALSE))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      isTRUE(try(rstudioapi::isAvailable(), silent = TRUE))) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (is.character(p) && length(p) == 1L && nzchar(p))
      return(normalizePath(dirname(p), winslash = "/", mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

find_data_file <- function(fname, here) {
  roots <- unique(c(here, file.path(here, "data"), dirname(here),
                    file.path(dirname(here), "data"),
                    getwd(), file.path(getwd(), "data"), dirname(getwd())))
  for (r in roots) {
    p <- file.path(r, fname)
    if (file.exists(p)) return(normalizePath(p, winslash = "/"))
  }
  loose <- list.files(unique(c(here, dirname(here))),
                      pattern = "^diltiazem.*\\.csv$", recursive = TRUE,
                      full.names = TRUE, ignore.case = TRUE)
  if (length(loose)) return(normalizePath(loose[1], winslash = "/"))
  stop("Data file '", fname, "' not found.\n",
       "  Folders searched: ", paste(roots, collapse = " | "), call. = FALSE)
}

HERE    <- script_dir()
fpath   <- find_data_file(DATA_NAME, HERE)
OUT_DIR <- file.path(HERE, "output_v2")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!nzchar(getOption("repos")[["CRAN"]]) ||
    identical(unname(getOption("repos")[["CRAN"]]), "@CRAN@"))
  options(repos = c(CRAN = "https://cloud.r-project.org"))

## nnet is needed because WeightIt fits the multinomial propensity model
## for a three-level treatment through nnet::multinom.
need <- c("WeightIt", "survey", "EValue", "nnet")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  message("Installing missing package(s): ", paste(miss, collapse = ", "))
  try(utils::install.packages(miss), silent = TRUE)   # no type = "binary"
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
}
if (length(miss)) stop("Missing package(s): ", paste(miss, collapse = ", "),
                       "\n  Install them with install.packages(), then re-run.",
                       call. = FALSE)
suppressPackageStartupMessages({ library(WeightIt); library(survey); library(EValue) })

LOG <- character(0)
say <- function(...) { txt <- paste0(...); LOG <<- c(LOG, txt)
                       cat(txt, "\n", sep = ""); invisible(NULL) }
rule  <- function() say(strrep("-", 74))
FAILS <- character(0)

say("E-VALUE LOG -- DILTIAZEM RVR-AF")
say("Run: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("R: ", R.version.string, " | ", Sys.info()[["sysname"]])
say("Scope: weighted risk ratios and E-values, primary outcome only")
say(strrep("=", 74))


## ---------------------------------------------------------------------
## 2. READ THE DATA
##    The archived file is semicolon-separated, uses the comma as the
##    decimal mark, and is pure ASCII apart from ONE byte: a Turkish
##    letter inside the unused column name "@yari_doz".
##
##    That byte matters. fileEncoding = "latin1" works on a Turkish
##    Windows machine but makes R abort the re-encoding connection at
##    that byte on Linux (verified on R 4.3.3: 19 columns, 0 rows). So
##    the file is read as raw bytes instead, any non-ASCII byte is
##    replaced by ".", and the text is parsed from memory. The same bytes
##    then give the same table on every operating system and locale.
##
##    check.names = FALSE keeps the column names exactly as they are in
##    the file; one of them contains a comma (LD_0,10_HR), which is why
##    the comma can never be used as the field separator here.
## ---------------------------------------------------------------------
say("Data file: ", fpath)

read_archive <- function(path) {
  con <- file(path, open = "rb"); on.exit(close(con), add = TRUE)
  bytes <- readBin(con, what = "raw", n = file.size(path) + 8L)
  bytes <- bytes[bytes != as.raw(0L)]
  if (length(bytes) >= 3L &&
      identical(bytes[1:3], as.raw(c(239L, 187L, 191L))))   # UTF-8 BOM
    bytes <- bytes[-(1:3)]
  bytes[bytes >= as.raw(128L)] <- charToRaw(".")
  lines <- strsplit(rawToChar(bytes), "\r\n|\n|\r")[[1]]
  utils::read.table(text = lines, header = TRUE, sep = ";", dec = ".",
                    quote = "\"", comment.char = "", check.names = FALSE,
                    fill = TRUE, stringsAsFactors = FALSE)
}
d <- read_archive(fpath)
if (ncol(d) < 40L || !all(c("Patient_NB", "Age", "Pulse_Control") %in% names(d)))
  stop("The data file did not parse as expected (", ncol(d), " columns, ",
       nrow(d), " rows).", call. = FALSE)
say("Read byte-wise, encoding-independent; separator ';', decimal comma per column")
say("Rows: ", nrow(d), "   Columns: ", ncol(d))
rule()


## ---------------------------------------------------------------------
## 3. HELPERS
## ---------------------------------------------------------------------
pick <- function(cands, what) {
  h <- cands[cands %in% names(d)]
  if (length(h) == 0L) stop("Column not found for '", what, "'.\n  Tried: ",
                            paste(cands, collapse = ", "), "\n  Available: ",
                            paste(names(d), collapse = ", "), call. = FALSE)
  h[1]
}
as_num <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  y <- gsub(" ", "", as.character(x), fixed = TRUE)
  y <- gsub(",", ".", y, fixed = TRUE)          # decimal comma -> decimal point
  suppressWarnings(as.numeric(y))
}
chr <- function(x) trimws(as.character(x))

check_int <- function(label, got, want) {
  ok <- isTRUE(all.equal(as.integer(got), as.integer(want)))
  if (!ok) FAILS <<- c(FAILS, label)
  say(sprintf("  %-26s %-22s expected %-22s %s", label,
              paste(got, collapse = " / "), paste(want, collapse = " / "),
              if (ok) "OK" else "MISMATCH"))
  ok
}
check_num <- function(label, got, want, tol = 0.02) {
  ok <- all(abs(got - want) < tol)
  if (!ok) FAILS <<- c(FAILS, label)
  say(sprintf("  %-26s %-22s expected %-22s %s", label,
              paste(sprintf("%.2f", got), collapse = " / "),
              paste(sprintf("%.2f", want), collapse = " / "),
              if (ok) "OK" else "MISMATCH"))
  ok
}


## ---------------------------------------------------------------------
## 4. TREATMENT GROUP
## ---------------------------------------------------------------------
V_GRP <- pick(c("treatment_group", "Treatment_Group"), "treatment group")
gv <- chr(d[[V_GRP]])
LV <- c("Low_Dose", "Standard_Bolus", "Standard_Infusion")

if (all(LV %in% unique(gv))) {
  d$GRP_F <- factor(gv, levels = LV)
  say("Group levels matched by name: ", paste(LV, collapse = ", "))
} else {
  uq <- sort(unique(gv[!is.na(gv) & gv != ""]))
  if (length(uq) != 3L) stop("Expected 3 treatment groups, found: ",
                             paste(uq, collapse = ", "), call. = FALSE)
  ord <- uq[order(-as.integer(table(gv)[uq]))]   ## 105 > 59 > 50, all distinct
  d$GRP_F <- factor(gv, levels = ord)
  LV <- ord
  say("Group levels matched by size (the file codes them numerically):")
  say("  ", paste(sprintf("%s -> n=%d", ord, as.integer(table(gv)[ord])),
                  collapse = " | "))
}


## ---------------------------------------------------------------------
## 5. BUILD THE NINE COVARIATES AND THE OUTCOME FROM THE RAW COLUMNS
## ---------------------------------------------------------------------
V_PULSE_C <- pick(c("Pulse_Control"), "10-minute heart rate")
V_VRR_REC <- pick(c("First_Dose_Ventricular_Rate_Response"), "recorded response")
V_GEN  <- pick(c("Gender"), "sex")
V_CAD  <- pick(c("CAD"), "coronary artery disease")
V_NAF  <- pick(c("New_Onset_AF"), "new-onset AF")
V_AA   <- pick(c("Use_of_Antiarrhytmic", "Use_of_Antiarrhythmic"), "antiarrhythmic use")
V_WAA  <- pick(c("Which_antiarrhytmic", "Which_antiarrhythmic"), "antiarrhythmic class")

for (v in c("Age", "Pulse_Admission", "Admission_SP", "Admission_MAP", "Weight", V_PULSE_C)) {
  if (!(v %in% names(d))) stop("Numeric column missing: ", v, call. = FALSE)
  d[[v]] <- as_num(d[[v]])
}

## Primary outcome: heart rate <= 100 bpm at 10 minutes, per the protocol
d$OUT_Y <- as.integer(d[[V_PULSE_C]] <= 100)

## Binary covariate: take whichever level reproduces the published counts.
## Which level is called "1" does not change the propensity model, only the
## sign of its coefficient, so the weights are the same either way.
mk_bin <- function(col, want, label) {
  x  <- chr(d[[col]])
  lv <- sort(unique(x[!is.na(x) & x != ""]))
  matched <- NA_character_
  for (l in lv) {
    cn <- as.integer(tapply(x == l, d$GRP_F, function(z) sum(z, na.rm = TRUE))[LV])
    if (isTRUE(all.equal(cn, as.integer(want)))) matched <- l
  }
  if (is.na(matched)) {
    FAILS <<- c(FAILS, label)
    say(sprintf("  %-26s levels {%s} -- NO level reproduces %s   MISMATCH",
                label, paste(lv, collapse = ","), paste(want, collapse = "/")))
    matched <- lv[length(lv)]
  } else {
    say(sprintf("  %-26s level '%s' = 1  -> %-18s expected %-18s OK",
                label, matched,
                paste(as.integer(tapply(x == matched, d$GRP_F,
                      function(z) sum(z, na.rm = TRUE))[LV]), collapse = " / "),
                paste(want, collapse = " / ")))
  }
  as.integer(x == matched)
}

say("DERIVED VARIABLES")
d$Gender_bin       <- mk_bin(V_GEN, EXP$male,  "Gender_bin (male)")
d$CAD_bin          <- mk_bin(V_CAD, EXP$cad,   "CAD_bin")
d$New_Onset_AF_bin <- mk_bin(V_NAF, EXP$newaf, "New_Onset_AF_bin")

## Chronic rate-control medication: try the yes/no flag first; if it does
## not reproduce 80/48/37, fall back to the drug-class column, where the
## "none" category is identified by its total of 49 patients (which does
## not depend on how the categories happen to be coded).
rc <- mk_bin(V_AA, EXP$ratectl, "RateControl_bin (flag)")
if ("RateControl_bin (flag)" %in% FAILS) {
  FAILS <- setdiff(FAILS, "RateControl_bin (flag)")
  w  <- chr(d[[V_WAA]])
  tb <- table(w)
  none_lv <- names(tb)[which(as.integer(tb) == (sum(EXP$n) - sum(EXP$ratectl)))]
  if (length(none_lv) == 1L) {
    rc <- as.integer(w != none_lv)
    cn <- as.integer(tapply(rc == 1, d$GRP_F, function(z) sum(z, na.rm = TRUE))[LV])
    ok <- isTRUE(all.equal(cn, as.integer(EXP$ratectl)))
    if (!ok) FAILS <- c(FAILS, "RateControl_bin (class)")
    say(sprintf("  %-26s 'none' level = '%s'  -> %-18s expected %-18s %s",
                "RateControl_bin (class)", none_lv, paste(cn, collapse = " / "),
                paste(EXP$ratectl, collapse = " / "), if (ok) "OK" else "MISMATCH"))
  } else {
    FAILS <- c(FAILS, "RateControl_bin")
    say("  RateControl_bin            could not be derived from either column   MISMATCH")
  }
}
d$RateControl_bin <- rc

COVS <- c("Age", "Gender_bin", "CAD_bin", "New_Onset_AF_bin", "RateControl_bin",
          "Pulse_Admission", "Admission_SP", "Admission_MAP", "Weight")

keep <- stats::complete.cases(d[, c(COVS, "OUT_Y")]) & !is.na(d$GRP_F)
if (sum(!keep) > 0L) say("  NOTE: ", sum(!keep), " row(s) dropped for missing values")
d <- d[keep, , drop = FALSE]
rule()


## ---------------------------------------------------------------------
## 6. VERIFY EVERY DERIVED QUANTITY AGAINST analysis_log_v2
## ---------------------------------------------------------------------
say("VERIFICATION vs analysis_log_v2")
gm <- function(v) as.numeric(tapply(d[[v]], d$GRP_F, mean, na.rm = TRUE)[LV])
check_int("n per group", as.integer(table(d$GRP_F)[LV]), EXP$n)
check_int("successes (HR <= 100)", as.integer(tapply(d$OUT_Y, d$GRP_F, sum)[LV]), EXP$success)
check_num("mean Age", gm("Age"), EXP$age)
check_num("mean Weight", gm("Weight"), EXP$weight)
check_num("mean Pulse_Admission", gm("Pulse_Admission"), EXP$pulse)
check_num("mean Admission_SP", gm("Admission_SP"), EXP$sbp)
check_num("mean Admission_MAP", gm("Admission_MAP"), EXP$map)

## independent check: does the recorded response column agree with the rule?
vr <- chr(d[[V_VRR_REC]])
agree <- NA
for (l in sort(unique(vr[!is.na(vr) & vr != ""]))) {
  if (isTRUE(all.equal(as.integer(vr == l), as.integer(d$OUT_Y)))) agree <- l
}
say(sprintf("  %-26s %s", "recorded response agrees",
            if (!is.na(agree)) paste0("yes, level '", agree, "' = success   OK")
            else "NO -- the recorded response does not match HR <= 100   MISMATCH"))
if (is.na(agree)) FAILS <- c(FAILS, "recorded response")
rule()


## ---------------------------------------------------------------------
## 7. STABILISED IPTW (ATE)
## ---------------------------------------------------------------------
frm <- stats::as.formula(paste("GRP_F ~", paste(COVS, collapse = " + ")))
set.seed(20261002)
W <- tryCatch(WeightIt::weightit(frm, data = d, method = "glm",
                                 estimand = "ATE", stabilize = TRUE),
              error = function(e) WeightIt::weightit(frm, data = d, method = "ps",
                                 estimand = "ATE", stabilize = TRUE))
d$WT_ <- as.numeric(W$weights)
ess <- as.numeric(tapply(d$WT_, d$GRP_F, function(w) sum(w)^2 / sum(w^2))[LV])

say("IPTW: stabilised weights, ATE, multinomial logistic propensity model")
say(sprintf("  weight range %.3f - %.3f   (expected 0.396 - 2.385)",
            min(d$WT_), max(d$WT_)))
check_num("effective sample sizes", ess, c(102.02, 50.14, 42.26), tol = 0.5)
rule()


## ---------------------------------------------------------------------
## 8. WEIGHTED RISK RATIOS AND E-VALUES
##    The first-named group goes in the numerator, as in the manuscript.
## ---------------------------------------------------------------------
ZC <- stats::qnorm(0.975)
rr_pair <- function(g1, g2) {
  s <- d[d$GRP_F %in% c(g1, g2), , drop = FALSE]
  s$CMP_ <- factor(as.character(s$GRP_F), levels = c(g2, g1))   ## g2 = reference
  des <- survey::svydesign(ids = ~1, weights = ~WT_, data = s)
  fit <- survey::svyglm(OUT_Y ~ CMP_, design = des,
                        family = stats::quasipoisson(link = "log"))
  cf <- summary(fit)$coefficients
  b <- cf[2, 1]; se <- cf[2, 2]
  list(est = exp(b), lo = exp(b - ZC * se), hi = exp(b + ZC * se), p = cf[2, 4])
}
ev_manual <- function(rr) {
  if (is.na(rr) || rr <= 0) return(NA_real_)
  r <- if (rr < 1) 1 / rr else rr
  r + sqrt(r * (r - 1))
}

PAIRS <- list(c(LV[1], LV[2]), c(LV[1], LV[3]), c(LV[2], LV[3]))
res <- NULL
say("WEIGHTED RISK RATIOS AND E-VALUES (primary outcome)")
say("")
for (pr in PAIRS) {
  lab <- paste(pr[1], "vs", pr[2])
  r   <- rr_pair(pr[1], pr[2])
  ev  <- EValue::evalues.RR(est = r$est, lo = r$lo, hi = r$hi)
  evp <- as.numeric(ev["E-values", "point"])
  evl <- as.numeric(ev["E-values", "lower"])
  if (is.na(evl)) evl <- as.numeric(ev["E-values", "upper"])
  crosses <- (r$lo <= 1 && r$hi >= 1)
  if (crosses || is.na(evl)) evl <- 1
  say(lab)
  say(sprintf("  RR            %.4f  (%.4f - %.4f)   p = %.4g", r$est, r$lo, r$hi, r$p))
  say(sprintf("  E-value       point %.4f   CI limit %.4f%s", evp, evl,
              if (crosses) "   [CI includes 1, so the limit E-value is 1]" else ""))
  say(sprintf("  manual check  point %.4f   (package minus manual = %.6f)",
              ev_manual(r$est), evp - ev_manual(r$est)))
  say(sprintf("  FOR THE TABLE  RR %.2f (%.2f-%.2f) | E %.2f | E(limit) %.2f",
              r$est, r$lo, r$hi, evp, evl))
  say("")
  res <- rbind(res, data.frame(comparison = lab, rr = r$est, lo = r$lo, hi = r$hi,
                               ev_point = evp, ev_limit = evl, stringsAsFactors = FALSE))
}
rule()


## ---------------------------------------------------------------------
## 9. CROSS-CHECK THE TWO PUBLISHED RISK RATIOS
## ---------------------------------------------------------------------
TOL <- 0.001
say("CROSS-CHECK vs analysis_log_v2")
rr_ok <- logical(0)
for (i in 1:2) {
  rf <- REF[[i]]
  de <- abs(res$rr[i] - rf["est"])
  rr_ok <- c(rr_ok, de < TOL)
  say("  ", res$comparison[i])
  say(sprintf("    RR  here %.7f | log %.7f | diff %.7f", res$rr[i], rf["est"], de))
  say(sprintf("    CI  here %.4f-%.4f | log %.4f-%.4f", res$lo[i], res$hi[i],
              rf["lo"], rf["hi"]))
}
say("")
if (length(FAILS) > 0L)
  say("  Verification items that did NOT reproduce: ", paste(unique(FAILS), collapse = ", "))
if (all(rr_ok) && length(FAILS) == 0L) {
  say("  VERDICT: PASS -- every checked quantity and both published risk ratios")
  say("  reproduce. The third comparison and its E-value come from the same")
  say("  model and can be reported in Table 5.")
} else if (all(rr_ok)) {
  say("  VERDICT: PASS WITH NOTES -- both published risk ratios reproduce, but")
  say("  the items listed above did not. Read those lines before using the")
  say("  third E-value.")
} else {
  say("  VERDICT: FAIL -- the published risk ratios are NOT reproduced.")
  say("  DO NOT report the third RR or its E-value from this run. Keep this log;")
  say("  the covariate construction or the propensity model differs.")
}
rule()


## ---------------------------------------------------------------------
## 10. SUMMARY BLOCK FOR TABLE 5
## ---------------------------------------------------------------------
say("TABLE 5 -- SENSITIVITY TO UNMEASURED CONFOUNDING")
say("")
say(sprintf("%-44s %-18s %-8s %-8s", "Comparison", "RR (95% CI)", "E point", "E limit"))
for (i in seq_len(nrow(res))) {
  say(sprintf("%-44s %-18s %-8s %-8s", res$comparison[i],
              sprintf("%.2f (%.2f-%.2f)", res$rr[i], res$lo[i], res$hi[i]),
              sprintf("%.2f", res$ev_point[i]), sprintf("%.2f", res$ev_limit[i])))
}
say("")
say("E-value method: VanderWeele & Ding 2017. The risk ratio was estimated")
say("directly; the sqrt(OR) approximation was not used.")
say(strrep("=", 74))
say("R: ", R.version.string)
for (p in c("WeightIt", "survey", "EValue", "nnet")) {
  v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed")
  say("  ", p, " ", v)
}
say("END OF LOG")


## ---------------------------------------------------------------------
## 11. WRITE THE LOG -- forced to plain ASCII
## ---------------------------------------------------------------------
out <- iconv(LOG, from = "", to = "ASCII", sub = "?")
out[is.na(out)] <- "?"
con <- file(file.path(OUT_DIR, LOGFILE), open = "wt")
writeLines(out, con, useBytes = TRUE)
close(con)
cat("\nWritten: ", file.path(OUT_DIR, LOGFILE), "\n", sep = "")
