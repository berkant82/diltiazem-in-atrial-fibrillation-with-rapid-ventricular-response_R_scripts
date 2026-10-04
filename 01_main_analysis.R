## =====================================================================
## DILTIAZEM FOR ATRIAL FIBRILLATION WITH RAPID VENTRICULAR RESPONSE
## Script 1 of 3  --  MAIN ANALYSIS
## =====================================================================
## This script reproduces every number reported in the manuscript in a
## single run: dose descriptives, baseline table, primary and secondary
## outcomes, infusion kinetics, regression within the low-dose group,
## propensity-score analyses, E-values and the sensitivity analyses.
##
## Primary outcome ........ heart rate <= 100 bpm, 10 min after treatment
## Primary adjusted model . stabilised IPTW (ATE); body WEIGHT is in the
##                          propensity model
## Secondary .............. 1:1 nearest-neighbour propensity-score
##                          matching (caliper 0.2 SD on the logit scale)
## Supplement ............. BMI model, entropy balancing, optimal full
##                          matching, Firth correction, balance diagnostics
##
## HOW TO RUN
##   Put this script in the same folder as (or in a sub-folder next to)
##   the data file "diltiazem_doz1set_eng_.csv" and run it. Nothing has
##   to be edited: the script finds its own folder, finds the data file,
##   creates an "output_v2" folder beside itself and writes everything
##   there. It behaves identically on Windows, macOS and Linux and does
##   not depend on the system language or the decimal separator of the
##   operating system.
##
##   From a terminal:   Rscript 01_main_analysis.R
##   From RStudio:      open the file and press "Source"
##   From an R console: source("01_main_analysis.R")
##
## OUTPUT
##   output_v2/analysis_log_v2.txt    full transcript of the run
##   output_v2/table2_first_dose.csv
##   output_v2/table2b_cumulative_dose.csv
##   output_v2/loveplot_*.png, output_v2/overlap_*.png
##
## Written for R >= 4.1. Tested with R 4.3.
## =====================================================================


## ---- 0. PORTABLE SETUP ----------------------------------------------
## Nothing below this block needs to be edited by the user.

options(stringsAsFactors = FALSE, OutDec = ".", warn = 1)

## Sorting and messages must not depend on the operating system language,
## otherwise factor levels could come out in a different order on a
## Turkish, German or Japanese machine and the tables would not match.
invisible(suppressWarnings(Sys.setlocale("LC_COLLATE", "C")))
Sys.setenv(LANGUAGE = "en")

## Close any output redirection left over from an interrupted earlier run,
## and make sure an error does not leave the log file open.
while (sink.number() > 0) sink()
options(error = function() {
  while (sink.number() > 0) sink()
  if (!interactive()) quit(save = "no", status = 1, runLast = FALSE)
})

## Locate the folder this script lives in, whichever way it was started
## (Rscript, source(), the RStudio "Source" button, or an R console).
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
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                  error = function(e) "")
    if (is.character(p) && length(p) == 1L && nzchar(p))
      return(normalizePath(dirname(p), winslash = "/", mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

## Look for the data file in the obvious places: beside the script, one
## level up (scripts in english/ or turkish/, data at the top), in a
## "data" sub-folder, and in the current working directory.
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
       "  Put it in the same folder as this script, or one level above it.\n",
       "  Folders searched: ", paste(roots, collapse = " | "), call. = FALSE)
}

HERE      <- script_dir()
DATA_FILE <- find_data_file("diltiazem_doz1set_eng_.csv", HERE)
OUT_DIR   <- file.path(HERE, "output_v2")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(2026)

## Install any missing packages from CRAN. Do NOT pass type = "binary":
## that argument is a hard error on Linux and is unnecessary elsewhere.
if (!nzchar(getOption("repos")[["CRAN"]]) ||
    identical(unname(getOption("repos")[["CRAN"]]), "@CRAN@"))
  options(repos = c(CRAN = "https://cloud.r-project.org"))

ensure_packages <- function(pkgs, required = TRUE) {
  need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(need)) {
    message("Installing missing package(s): ", paste(need, collapse = ", "))
    try(utils::install.packages(need), silent = TRUE)
  }
  ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  if (required && !all(ok))
    stop("These packages are required but could not be installed: ",
         paste(pkgs[!ok], collapse = ", "),
         "\n  Install them manually with install.packages(), then re-run.",
         call. = FALSE)
  invisible(lapply(pkgs[ok], function(p)
    suppressPackageStartupMessages(library(p, character.only = TRUE))))
  ok
}

## nnet is listed because WeightIt fits the multinomial propensity model
## for a three-level treatment through nnet::multinom.
ensure_packages(c("MatchIt", "WeightIt", "cobalt", "survey", "tableone", "nnet"))
has_logistf <- all(ensure_packages("logistf", required = FALSE))
has_evalue  <- all(ensure_packages("EValue",  required = FALSE))

## Save a plot without assuming that a screen graphics device exists
## (a bare Linux server has no X11, so png() alone can fail there).
save_plot <- function(p, file, width = 9, height = 6, dpi = 100) {
  if (requireNamespace("ggplot2", quietly = TRUE) && inherits(p, "ggplot")) {
    ok <- tryCatch({
      ggplot2::ggsave(file, p, width = width, height = height, dpi = dpi)
      TRUE }, error = function(e) FALSE)
    if (ok) return(invisible(TRUE))
  }
  ok <- tryCatch({
    grDevices::png(file, width = width * dpi, height = height * dpi, res = dpi)
    print(p)
    grDevices::dev.off()
    TRUE }, error = function(e) { try(grDevices::dev.off(), silent = TRUE); FALSE })
  if (!ok) message("NOTE: no usable graphics device, skipped ", basename(file))
  invisible(ok)
}

sink(file.path(OUT_DIR, "analysis_log_v2.txt"), split = TRUE)
cat("DILTIAZEM RVR-AF - ANALYSIS LOG v2\n")
cat("Run on:", format(Sys.time()), "\n")
cat("R:", R.version.string, "|", Sys.info()[["sysname"]], "\n")
cat("Script folder:", HERE, "\n")
cat("Data file    :", DATA_FILE, "\n")
cat("Output folder:", OUT_DIR, "\n")
cat("Primary outcome: heart rate <= 100 bpm, 10 min after start of treatment\n")
cat("Primary adjusted analysis: stabilised IPTW (ATE), body weight in the PS model\n")
cat(strrep("=", 72), "\n")


## ---- 1. DATA AND INTEGRITY CHECKS -----------------------------------
## The archived file is semicolon-separated, uses the comma as the decimal
## mark, and is pure ASCII apart from ONE byte: a Turkish letter inside the
## column name "@yari_doz", which no analysis uses.
##
## That single byte matters. Passing fileEncoding = "latin1" to read.csv()
## works on a Turkish Windows machine but makes R abort the re-encoding
## connection at that byte on Linux (verified on R 4.3.3: the file comes
## back as 19 columns and 0 rows). So the file is not read through an
## encoding connection at all. It is read as raw bytes, any byte outside
## ASCII is replaced by "." - which is exactly what make.names() would do
## to it anyway - and the result is parsed from text. Identical bytes then
## give an identical data frame on every operating system and in every
## locale, with nothing for the user to configure.
read_archive_csv <- function(path, check_names = TRUE) {
  con <- file(path, open = "rb"); on.exit(close(con), add = TRUE)
  bytes <- readBin(con, what = "raw", n = file.size(path) + 8L)
  bytes <- bytes[bytes != as.raw(0L)]
  if (length(bytes) >= 3L &&
      identical(bytes[1:3], as.raw(c(239L, 187L, 191L))))   # UTF-8 BOM
    bytes <- bytes[-(1:3)]
  bytes[bytes >= as.raw(128L)] <- charToRaw(".")
  lines <- strsplit(rawToChar(bytes), "\r\n|\n|\r")[[1]]
  utils::read.csv(text = lines, sep = ";", header = TRUE,
                  check.names = check_names, stringsAsFactors = FALSE)
}
df <- read_archive_csv(DATA_FILE)
cat("File read byte-wise (encoding-independent) | rows:", nrow(df),
    "| columns:", ncol(df), "\n")
if (ncol(df) < 40L || !all(c("Patient_NB", "Age", "Pulse_Control") %in% names(df)))
  stop("The data file did not parse as expected (", ncol(df), " columns, ",
       nrow(df), " rows). Is it the archived semicolon-separated file?",
       call. = FALSE)

names(df)[4] <- "treatment_group"

## Normalise the column names. R turns a leading "@" into "." or into "X."
## depending on the version, so force a single spelling. The names are
## already ASCII after the byte-wise read; the iconv() call is a safety
## net in case someone substitutes a differently encoded file.
nm <- iconv(names(df), from = "", to = "ASCII", sub = ".")
nm[is.na(nm)] <- names(df)[is.na(nm)]
names(df) <- sub("^[.]", "X.", make.names(nm, unique = TRUE))

## Does the file really contain the columns the code uses?
needed <- c("Patient_NB","Age","Gender","treatment_group","CAD","New_Onset_AF",
            "Use_of_Antiarrhytmic","Which_antiarrhytmic","Pulse_Admission",
            "Pulse_Control","First_Dose_Ventricular_Rate_Response","Result",
            "Admission_SP","Admission_MAP","Control_Systolic_Blood_Pressure",
            "Control_MAP","Low_dose","aditional_dose","additional_medication",
            "Secound_Dose_Ventricular_Rate_Response",
            "X.Beta_Bloker_Ventricular_Rate_Response","X.0.35_mg_dose",
            "Complication","Weight","Height","BMI","Ordered_Dose",
            "X.150sec_VRR","X.300sec_VRR","X.450sec_VRR","X.600sec_VRR",
            "X.150sec_Pulse","X.300sec_Pulse","X.450sec_Pulse","X.600sec_Pulse")
absent <- setdiff(needed, names(df))
if (length(absent)) {
  cat("!!! MISSING COLUMNS:\n"); print(absent)
  cat("\nNames found in the file:\n"); print(names(df))
  stop("Column names do not match - see the list above.", call. = FALSE)
}
cat("Column-name check OK\n")

## Group codes in the file: 3 = low dose, 2 = standard bolus,
## 1 = standard infusion.
df$treatment_group <- factor(df$treatment_group, levels = c(3, 2, 1),
                             labels = c("Low_Dose", "Standard_Bolus", "Standard_Infusion"))

## Decimal comma -> decimal point. These columns are stored as text in the
## archive because they were exported from a Turkish-locale spreadsheet.
numeric_cols <- c("Admission_MAP","Control_MAP","BMI","Weight","Height","Ordered_Dose",
                  "After_First_Dose_HR","After_Sec_Dose_Hr","After_BB_HR")
for (v in numeric_cols) if (v %in% names(df))
  df[[v]] <- as.numeric(gsub(",", ".", as.character(df[[v]])))

## SUCCESS THRESHOLD: heart rate <= 100 bpm
THRESHOLD <- 100
df$VRR_binary   <- ifelse(df$First_Dose_Ventricular_Rate_Response == 1, 1, 0)
df$Comp_binary  <- ifelse(df$Complication == 1, 1, 0)
df$Adm_binary   <- ifelse(df$Result == 2, 1, 0)
df$Gender_bin   <- ifelse(df$Gender == 1, 1, 0)
df$CAD_bin      <- ifelse(df$CAD == 1, 1, 0)
df$New_Onset_AF_bin  <- ifelse(df$New_Onset_AF == 1, 1, 0)
df$RateControl_bin   <- ifelse(df$Use_of_Antiarrhytmic == 1, 1, 0)
df$RateControl_class <- factor(df$Which_antiarrhytmic, levels = c(1, 2, 4, 3),
                               labels = c("CCB", "Beta_blocker", "Digoxin", "None"))

## --- integrity: does the recorded response flag agree with the rule? ---
agrees <- (df$Pulse_Control <= THRESHOLD) == (df$VRR_binary == 1)
if (!all(agrees, na.rm = TRUE)) {
  cat("\n!!! WARNING: patients where the recorded response disagrees with HR<=",
      THRESHOLD, ":\n", sep = "")
  print(df[!agrees, c("Patient_NB","treatment_group","Pulse_Control",
                      "First_Dose_Ventricular_Rate_Response")])
} else cat("\nINTEGRITY OK: the recorded response matches HR<=", THRESHOLD,
           " for every patient\n", sep = "")

## --- integrity: cohort size and group sizes ---
observed_n <- as.integer(table(df$treatment_group))
if (nrow(df) != 214L || !identical(observed_n, c(105L, 59L, 50L)))
  stop("Cohort does not match the published one. Expected n = 214 ",
       "(105 / 59 / 50), found n = ", nrow(df), " (",
       paste(observed_n, collapse = " / "), ").\n",
       "  The data file is probably not the archived version.", call. = FALSE)
cat("n = 214 (105 / 59 / 50)\n")
cat("Successes:", paste(sprintf("%s %d", levels(df$treatment_group),
    tapply(df$VRR_binary, df$treatment_group, sum)), collapse = " | "), "\n")

analysis_vars <- c("Age","Gender_bin","CAD_bin","New_Onset_AF_bin","RateControl_bin",
                   "Weight","BMI","Pulse_Admission","Admission_SP","Admission_MAP",
                   "Pulse_Control","Control_Systolic_Blood_Pressure","Control_MAP",
                   "Ordered_Dose","VRR_binary","Comp_binary","Adm_binary")
n_missing <- sapply(df[analysis_vars], function(x) sum(is.na(x)))
cat("Missing data:", if (all(n_missing == 0)) "none\n" else
    paste(names(n_missing[n_missing > 0]), n_missing[n_missing > 0], collapse = "; "), "\n")


## ---- 2. HELPER FUNCTIONS --------------------------------------------
group_pairs <- list(c("Low_Dose","Standard_Bolus"),
                    c("Low_Dose","Standard_Infusion"),
                    c("Standard_Bolus","Standard_Infusion"))
PS_COVARIATES <- c("Age","Gender_bin","CAD_bin","New_Onset_AF_bin","RateControl_bin",
                   "Pulse_Admission","Admission_SP","Admission_MAP")

## OR = odds(g1) / odds(g2), with a Wald confidence interval.
or_wald <- function(data, g1, g2, outcome) {
  d <- droplevels(subset(data, treatment_group %in% c(g1, g2)))
  d$tr <- ifelse(d$treatment_group == g1, 1, 0)
  co <- summary(glm(reformulate("tr", outcome), data = d,
                    family = binomial))$coefficients["tr", ]
  sprintf("%-38s OR %.2f (%.2f-%.2f), p=%.4g", paste(g1, "vs", g2),
          exp(co[1]), exp(co[1]-1.96*co[2]), exp(co[1]+1.96*co[2]), co[4])
}

## Odds ratio on the matched sample: conditional maximum likelihood with
## an exact confidence interval (fisher.test). Note that this treats the
## matched observations as two independent groups; the pair-conditional
## estimate is produced separately by script 02.
or_exact <- function(data, g1, g2, outcome) {
  m <- rbind(c(sum(data[[outcome]][data$treatment_group == g1] == 1),
               sum(data[[outcome]][data$treatment_group == g1] == 0)),
             c(sum(data[[outcome]][data$treatment_group == g2] == 1),
               sum(data[[outcome]][data$treatment_group == g2] == 0)))
  if (any(m == 0)) return(list(zero = TRUE, txt = sprintf(
    "%d/%d vs %d/%d - NOT ESTIMABLE (zero cell)",
    m[1,1], sum(m[1,]), m[2,1], sum(m[2,]))))
  ft <- fisher.test(m)
  list(zero = FALSE, txt = sprintf("OR %.2f (%.2f-%.2f), p=%.4g",
       ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))
}

## Welch's heteroscedastic one-way ANOVA, written out by hand because
## oneway.test() raises a formula-checking error in some R versions.
##   lambda = 3 * S / (k^2 - 1),  where S = sum[(1 - w/W)^2 / (n - 1)]
##   F      = A / (1 + 2(k-2)S/(k^2-1)) = A / (1 + 2(k-2)*lambda/3)
##   df2    = (k^2 - 1) / (3S) = 1 / lambda
welch <- function(y, g) {
  ok <- !is.na(y) & !is.na(g); y <- y[ok]; g <- droplevels(factor(g[ok]))
  n <- tapply(y,g,length); m <- tapply(y,g,mean); v <- tapply(y,g,var); k <- nlevels(g)
  w <- n/v; W <- sum(w); mw <- sum(w*m)/W
  A <- sum(w*(m-mw)^2)/(k-1); lam <- 3*sum((1-w/W)^2/(n-1))/(k^2-1)
  Fs <- A/(1+2*(k-2)*lam/3); df2 <- 1/lam
  list(F = Fs, df1 = k-1, df2 = df2, p = pf(Fs, k-1, df2, lower.tail = FALSE))
}

## Weighted odds ratio from a survey design (robust standard errors).
or_weighted <- function(design, g1, g2, outcome) {
  d <- subset(design, treatment_group %in% c(g1, g2))
  d <- update(d, tr = as.numeric(treatment_group == g1))
  co <- summary(svyglm(reformulate("tr", outcome), design = d,
                       family = quasibinomial))$coefficients["tr", ]
  sprintf("%-38s OR %.2f (%.2f-%.2f), p=%.4g", paste(g1, "vs", g2),
          exp(co[1]), exp(co[1]-1.96*co[2]), exp(co[1]+1.96*co[2]), co[4])
}

## Weighted risk ratio (log-link Poisson, robust SE). The E-value needs a
## risk ratio, not an odds ratio.
rr_weighted <- function(design, g1, g2, outcome) {
  d <- subset(design, treatment_group %in% c(g1, g2))
  d <- update(d, tr = as.numeric(treatment_group == g1))
  m <- svyglm(reformulate("tr", outcome), design = d,
              family = quasipoisson(link = "log"))
  co <- summary(m)$coefficients["tr", ]
  c(RR = unname(exp(co[1])), lo = unname(exp(co[1]-1.96*co[2])),
    hi = unname(exp(co[1]+1.96*co[2])), p = unname(co[4]))
}

## E-value (VanderWeele & Ding 2017), computed straight from the risk ratio.
e_value <- function(RR, lo, hi) {
  r <- unname(RR); L <- unname(lo); U <- unname(hi)
  if (r < 1) { nl <- 1/U; U <- 1/L; L <- nl; r <- 1/r }
  ev <- function(x) if (x <= 1) 1 else x + sqrt(x*(x-1))
  c(RR_direction = r, E_point = ev(r), E_CI = ev(L))
}

## Crude risk difference with a Wald confidence interval.
risk_difference <- function(data, g1, g2, outcome) {
  a <- data[[outcome]][data$treatment_group == g1]
  b <- data[[outcome]][data$treatment_group == g2]
  p1 <- mean(a); p2 <- mean(b); n1 <- length(a); n2 <- length(b)
  se <- sqrt(p1*(1-p1)/n1 + p2*(1-p2)/n2); d <- p1 - p2
  sprintf("%-38s RD %+.1f%% (%.1f to %.1f)", paste(g1, "vs", g2),
          100*d, 100*(d-1.96*se), 100*(d+1.96*se))
}


## ---- 3. DOSE ANALYSES -----------------------------------------------
cat("\n", strrep("=",72), "\n### DOSE ANALYSES\n", sep="")

df$protocol_mgkg <- ifelse(df$treatment_group == "Low_Dose", 0.15, 0.25)
cat("\n-- Ordered_Dose check (body weight x protocol coefficient) --\n")
deviation <- abs(df$Ordered_Dose - df$Weight * df$protocol_mgkg)
cat(sprintf("Exact agreement: %d/%d | largest deviation %.4f mg\n",
            sum(deviation < 0.001), nrow(df), max(deviation)))

cat("\n-- TABLE 2: first dose --\n")
dose_tab <- do.call(rbind, lapply(split(df, df$treatment_group), function(d) data.frame(
  n = nrow(d),
  weight = sprintf("%.1f +/- %.1f (%.0f-%.0f)", mean(d$Weight), sd(d$Weight),
                   min(d$Weight), max(d$Weight)),
  dose_mg = sprintf("%.2f +/- %.2f (%.2f-%.2f)", mean(d$Ordered_Dose), sd(d$Ordered_Dose),
                    min(d$Ordered_Dose), max(d$Ordered_Dose)),
  mg_kg = sprintf("%.2f (protocol)", unique(d$protocol_mgkg)))))
print(dose_tab)
write.csv(dose_tab, file.path(OUT_DIR, "table2_first_dose.csv"))

## rescue doses and cumulative exposure
low_dose_col <- trimws(as.character(df$Low_dose))
dose035      <- trimws(as.character(df$X.0.35_mg_dose))
df$extra_mgkg <- 0
df$extra_mgkg <- df$extra_mgkg +
  ifelse(df$treatment_group == "Low_Dose" & low_dose_col %in% c("2","3"), 0.10, 0)
df$extra_mgkg <- df$extra_mgkg + ifelse(dose035 %in% c("1","2"), 0.35, 0)
df$total_mgkg <- df$protocol_mgkg + df$extra_mgkg
df$total_mg   <- df$total_mgkg * df$Weight

cat("\n-- Rescue treatment --\n")
print(table(df$treatment_group, trimws(as.character(df$aditional_dose)),
            dnn = c("group","extra dose (1=yes)")))
cat("\n0.35 mg/kg rescue:\n")
print(table(df$treatment_group, dose035, dnn = c("group","0.35")))
cat("\nBeta blocker:\n")
print(table(df$treatment_group,
            trimws(as.character(df$X.Beta_Bloker_Ventricular_Rate_Response)),
            dnn = c("group","BB")))

cat("\n-- CUMULATIVE EXPOSURE (first dose + rescue) --\n")
cum_tab <- do.call(rbind, lapply(split(df, df$treatment_group), function(d) data.frame(
  n = nrow(d),
  first_mgkg = sprintf("%.2f", unique(d$protocol_mgkg)),
  total_mgkg = sprintf("%.3f +/- %.3f", mean(d$total_mgkg), sd(d$total_mgkg)),
  total_mg   = sprintf("%.1f +/- %.1f", mean(d$total_mg), sd(d$total_mg)),
  received_extra = sum(d$extra_mgkg > 0))))
print(cum_tab)
cat("\nCumulative mg/kg across groups (Welch):\n")
wk <- welch(df$total_mgkg, df$treatment_group)
cat(sprintf("  F=%.3f df=(%d, %.1f) p=%.4g\n", wk$F, wk$df1, wk$df2, wk$p))
write.csv(cum_tab, file.path(OUT_DIR, "table2b_cumulative_dose.csv"))


## ---- 4. TABLE 1: BASELINE + STANDARDISED DIFFERENCES -----------------
cat("\n", strrep("=",72), "\n### TABLE 1 - BASELINE\n", sep="")
t1 <- CreateTableOne(
  vars = c("Age","Gender_bin","CAD_bin","New_Onset_AF_bin","RateControl_bin",
           "Weight","BMI","Pulse_Admission","Admission_SP","Admission_MAP"),
  strata = "treatment_group", data = df,
  factorVars = c("Gender_bin","CAD_bin","New_Onset_AF_bin","RateControl_bin"))
## No p-value column: baseline p-values are not informative here and were
## removed at the reviewers' request (R2-M5). Balance is judged by SMD.
print(t1, smd = FALSE, test = FALSE)

cat("\n-- SMD before adjustment (largest of the pairwise comparisons, binary='std') --\n")
print(bal.tab(treatment_group ~ Age + Gender_bin + CAD_bin + New_Onset_AF_bin +
                RateControl_bin + Weight + BMI + Pulse_Admission +
                Admission_SP + Admission_MAP,
              data = df, un = TRUE, binary = "std", thresholds = c(m = .1)))
cat("\n-- Class of chronic rate-control drug --\n")
print(table(df$treatment_group, df$RateControl_class))
cat(sprintf("\nBMI range: %.1f - %.1f kg/m2 | BMI >= 30: %d patients\n",
            min(df$BMI), max(df$BMI), sum(df$BMI >= 30)))


## ---- 5. TABLE 3: PRIMARY OUTCOME -------------------------------------
cat("\n", strrep("=",72), "\n### TABLE 3 - PRIMARY OUTCOME (HR <= 100 bpm at 10 min)\n", sep="")
vrr <- table(df$treatment_group, df$VRR_binary)[, c("1","0")]
print(vrr)
cat(sprintf("Success: %s\n", paste(sprintf("%s %.1f%%", rownames(vrr),
      100*vrr[,1]/rowSums(vrr)), collapse=" | ")))
omnibus <- chisq.test(vrr)
cat(sprintf("Omnibus chi-square: X2=%.2f, df=%d, p=%.3g\n",
            omnibus$statistic, omnibus$parameter, omnibus$p.value))
cat("\nPairwise comparisons, Bonferroni corrected (chi-square with Yates correction):\n")
print(pairwise.prop.test(vrr[,1], rowSums(vrr), p.adjust.method = "bonferroni"))
cat("\n-- Unadjusted OR (Wald) --\n")
for (cf in group_pairs) cat(or_wald(df, cf[1], cf[2], "VRR_binary"), "\n")
cat("\n-- Crude risk difference --\n")
for (cf in group_pairs) cat(risk_difference(df, cf[1], cf[2], "VRR_binary"), "\n")


## ---- 6. TABLE 4: SECONDARY OUTCOMES ----------------------------------
cat("\n", strrep("=",72), "\n### TABLE 4 - SECONDARY OUTCOMES\n", sep="")
for (v in c("Pulse_Control","Control_Systolic_Blood_Pressure","Control_MAP")) {
  summ <- tapply(df[[v]], df$treatment_group, function(x)
    sprintf("%.1f +/- %.1f", mean(x), sd(x)))
  wa <- welch(df[[v]], df$treatment_group)
  bt <- bartlett.test(df[[v]], df$treatment_group)$p.value
  cat(sprintf("%-33s %s\n%-33s Welch F=%.3f df=(%d,%.1f) p=%.4g | Bartlett p=%.4g\n",
              v, paste(summ, collapse=" | "), "", wa$F, wa$df1, wa$df2, wa$p, bt))
}
cat("\n-- Adverse events --\n")
comp <- table(df$treatment_group, df$Comp_binary)[, c("1","0")]
print(comp)
cat("Omnibus Fisher p =", round(fisher.test(comp)$p.value, 4), "\n")
for (cf in group_pairs) cat(or_wald(df, cf[1], cf[2], "Comp_binary"), "\n")

cat("\n-- Breakdown by event type --\n")
events <- subset(df, Comp_binary == 1,
                 select = c(Patient_NB, treatment_group, Pulse_Control,
                            Control_Systolic_Blood_Pressure))
events$hypotension <- events$Control_Systolic_Blood_Pressure < 90
events$bradycardia <- events$Pulse_Control < 50
print(events, row.names = FALSE)
cat(sprintf("Total %d events | hypotension %d | bradycardia %d | other %d\n",
            nrow(events), sum(events$hypotension), sum(events$bradycardia),
            sum(!events$hypotension & !events$bradycardia)))

cat("\n-- Hospital admission (post hoc) --\n")
adm <- table(df$treatment_group, df$Adm_binary)[, c("1","0")]
print(adm); cat("Fisher p =", round(fisher.test(adm)$p.value, 4), "\n")


## ---- 7. INFUSION KINETICS --------------------------------------------
cat("\n", strrep("=",72), "\n### INFUSION KINETICS (standard infusion, n=50)\n", sep="")
si <- subset(df, treatment_group == "Standard_Infusion")
vrr_cols   <- c("X.150sec_VRR","X.300sec_VRR","X.450sec_VRR","X.600sec_VRR")
pulse_cols <- c("X.150sec_Pulse","X.300sec_Pulse","X.450sec_Pulse","X.600sec_Pulse")
for (k in pulse_cols) si[[k]] <- as.numeric(gsub(",", ".", as.character(si[[k]])))
cat("Source: the recorded response columns (not recomputed from heart rate)\n\n")
for (i in seq_along(vrr_cols)) {
  v  <- as.numeric(si[[vrr_cols[i]]])
  nb <- sum(si[[pulse_cols[i]]] <= THRESHOLD, na.rm = TRUE)
  cat(sprintf("%4d s: response column %2d/50 = %2.0f%% | from heart rate HR<=%d: %2d %s\n",
              i*150, sum(v == 1), 100*mean(v == 1), THRESHOLD, nb,
              ifelse(sum(v == 1) == nb, "", "<<< DIFFERENT")))
}
cat("\nPaired t-test against baseline, Bonferroni (4 comparisons):\n")
kp <- sapply(pulse_cols, function(k)
  t.test(si$Pulse_Admission, si[[k]], paired = TRUE)$p.value)
for (i in seq_along(kp))
  cat(sprintf("  %4d s: p=%.3g (Bonferroni %.3g)\n", i*150, kp[i],
              p.adjust(kp, "bonferroni")[i]))


## ---- 8. TABLE 5: REGRESSION WITHIN THE LOW-DOSE GROUP ----------------
cat("\n", strrep("=",72), "\n### TABLE 5 - PREDICTORS WITHIN THE LOW-DOSE GROUP\n", sep="")
ld <- subset(df, treatment_group == "Low_Dose")
m5 <- glm(VRR_binary ~ Age + Gender_bin + CAD_bin + New_Onset_AF_bin + RateControl_bin,
          data = ld, family = binomial)
co <- summary(m5)$coefficients
print(round(data.frame(OR = exp(co[,1]), lo = exp(co[,1]-1.96*co[,2]),
                       hi = exp(co[,1]+1.96*co[,2]), p = co[,4])[-1, ], 3))
cat(sprintf("n=%d, successes=%d, failures=%d (EPV = %.1f)\n",
            nrow(ld), sum(ld$VRR_binary), sum(1-ld$VRR_binary),
            sum(1-ld$VRR_binary)/5))
m0 <- glm(VRR_binary ~ 1, data = ld, family = binomial)
lr <- anova(m0, m5, test = "Chisq")
cat(sprintf("Model fit: LR chi-square = %.2f, df = %d, p = %.3f\n",
            lr$Deviance[2], lr$Df[2], lr$`Pr(>Chi)`[2]))
if (has_logistf) {
  f5 <- logistf(VRR_binary ~ Age + Gender_bin + CAD_bin + New_Onset_AF_bin +
                  RateControl_bin, data = ld)
  cat("Firth sensitivity analysis:\n")
  print(round(exp(cbind(OR = coef(f5), lo = f5$ci.lower, hi = f5$ci.upper))[-1, ], 3))
} else cat("NOTE: package 'logistf' unavailable, Firth analysis skipped.\n")


## ---- 9. TABLE 6: PROPENSITY-SCORE MATCHING AND IPTW -------------------
cat("\n", strrep("=",72), "\n### TABLE 6 - PROPENSITY SCORE ANALYSES\n", sep="")
ps_f <- as.formula(paste("tr ~", paste(c(PS_COVARIATES, "Weight"), collapse = " + ")))

cat("\n--- 9a. Matching (1:1 nearest neighbour, caliper 0.2 SD on the logit scale, ATT) ---\n")
cat("    In each pair the SECOND group named is treated as the 'treated'\n")
cat("    group, so the estimand is the ATT with respect to that group.\n")
cat("    Odds ratios are reported with the FIRST group in the numerator.\n")
for (cf in group_pairs) {
  d <- droplevels(subset(df, treatment_group %in% cf))
  d$tr <- ifelse(d$treatment_group == cf[2], 1, 0)
  ps <- matchit(ps_f, data = d, method = "nearest", distance = "glm",
                link = "linear.logit", caliper = 0.2, replace = FALSE)
  md <- match.data(ps)
  cat("\n", paste(cf[1], "vs", cf[2]), "- matched pairs:", sum(md$tr == 1), "\n")
  print(bal.tab(ps, un = TRUE, binary = "std", thresholds = c(m = .1)))
  bt <- bal.tab(ps, binary = "std", thresholds = c(m = .1))$Balance
  bt <- bt[rownames(bt) != "distance", ]
  cat(sprintf("  >> Covariates above the 0.1 threshold: %d / %d\n",
              sum(abs(bt$Diff.Adj) > 0.1, na.rm = TRUE), nrow(bt)))
  for (oc in c("VRR_binary","Comp_binary")) {
    r <- or_exact(md, cf[1], cf[2], oc)
    cat(sprintf("  %-12s %s\n", oc, r$txt))
    if (r$zero && has_logistf) {
      md$tr2 <- ifelse(md$treatment_group == cf[1], 1, 0)
      fm <- logistf(reformulate("tr2", oc), data = md)
      cat(sprintf("               Firth OR %.2f (%.2f-%.2f)\n", exp(coef(fm)["tr2"]),
                  exp(fm$ci.lower["tr2"]), exp(fm$ci.upper["tr2"])))
    }
  }
  save_plot(love.plot(ps, thresholds = .1, binary = "std",
                      title = paste("Covariate balance:", cf[1], "vs", cf[2])),
            file.path(OUT_DIR, paste0("loveplot_PSM_", cf[1], "_vs_", cf[2], ".png")))
  save_plot(bal.plot(ps, var.name = "distance", which = "both",
                     title = paste("PS overlap:", cf[1], "vs", cf[2])),
            file.path(OUT_DIR, paste0("overlap_", cf[1], "_vs_", cf[2], ".png")))
}

cat("\n--- 9b. IPTW (stabilised, ATE) - PRIMARY ADJUSTED ANALYSIS ---\n")
iptw <- weightit(as.formula(paste("treatment_group ~",
                 paste(c(PS_COVARIATES,"Weight"), collapse=" + "))),
                 data = df, method = "glm", estimand = "ATE", stabilize = TRUE)
print(summary(iptw))
print(bal.tab(iptw, un = TRUE, binary = "std", thresholds = c(m = .1)))
save_plot(love.plot(iptw, thresholds = .1, binary = "std",
                    title = "Covariate balance: IPTW (ATE, weight model)"),
          file.path(OUT_DIR, "loveplot_IPTW_primary.png"))

df$w <- iptw$weights
dsn  <- svydesign(ids = ~1, weights = ~w, data = df)
cat("\nIPTW odds ratios (robust SE):\n")
for (oc in c("VRR_binary","Comp_binary")) {
  cat(" --", oc, "--\n")
  for (cf in group_pairs) cat("  ", or_weighted(dsn, cf[1], cf[2], oc), "\n")
}
cat("\nIPTW risk ratios (log-Poisson, robust SE) - primary outcome:\n")
for (cf in group_pairs) {
  r <- rr_weighted(dsn, cf[1], cf[2], "VRR_binary")
  cat(sprintf("  %-38s RR %.2f (%.2f-%.2f), p=%.4g\n",
              paste(cf[1],"vs",cf[2]), r["RR"], r["lo"], r["hi"], r["p"]))
}


## ---- 10. E-VALUE ------------------------------------------------------
cat("\n", strrep("=",72), "\n### E-VALUE (from the weighted risk ratio)\n", sep="")
cat("Method: VanderWeele & Ding 2017. The risk ratio was estimated directly;\n")
cat("the sqrt(OR) approximation was NOT used (it is invalid for a common outcome).\n\n")
for (cf in group_pairs[1:2]) {
  r <- rr_weighted(dsn, cf[1], cf[2], "VRR_binary")
  e <- e_value(r["RR"], r["lo"], r["hi"])
  cat(sprintf("%s vs %s:\n  RR %.3f (%.3f-%.3f) -> E-value %.2f (at the CI limit %.2f)\n",
              cf[1], cf[2], r["RR"], r["lo"], r["hi"], e["E_point"], e["E_CI"]))
  if (has_evalue) {
    cat("  Cross-check with the EValue package:\n")
    print(evalues.RR(est = unname(r["RR"]), lo = unname(r["lo"]), hi = unname(r["hi"])))
  }
}


## ---- 11. SUPPLEMENT: SENSITIVITY ANALYSES -----------------------------
cat("\n", strrep("=",72), "\n### SUPPLEMENT - SENSITIVITY ANALYSES\n", sep="")

cat("\n--- S1. IPTW with BMI in the PS model instead of weight ---\n")
ib <- weightit(as.formula(paste("treatment_group ~",
               paste(c(PS_COVARIATES,"BMI"), collapse=" + "))),
               data = df, method = "glm", estimand = "ATE", stabilize = TRUE)
print(bal.tab(ib, binary = "std", thresholds = c(m = .1)))
df$wb <- ib$weights; dsb <- svydesign(ids = ~1, weights = ~wb, data = df)
for (oc in c("VRR_binary","Comp_binary")) { cat(" --", oc, "--\n")
  for (cf in group_pairs) cat("  ", or_weighted(dsb, cf[1], cf[2], oc), "\n") }

cat("\n--- S2. Entropy balancing (weight) ---\n")
ie <- weightit(as.formula(paste("treatment_group ~",
               paste(c(PS_COVARIATES,"Weight"), collapse=" + "))),
               data = df, method = "ebal", estimand = "ATE")
print(summary(ie)); print(bal.tab(ie, binary = "std", thresholds = c(m = .1)))
df$we <- ie$weights; dse <- svydesign(ids = ~1, weights = ~we, data = df)
for (oc in c("VRR_binary","Comp_binary")) { cat(" --", oc, "--\n")
  for (cf in group_pairs) cat("  ", or_weighted(dse, cf[1], cf[2], oc), "\n") }

cat("\n--- S3. Optimal full matching ---\n")
for (cf in group_pairs) {
  d <- droplevels(subset(df, treatment_group %in% cf))
  d$tr <- ifelse(d$treatment_group == cf[2], 1, 0)
  mf <- matchit(ps_f, data = d, method = "full", distance = "glm", estimand = "ATT")
  cat("\n", paste(cf[1],"vs",cf[2]), "\n")
  print(bal.tab(mf, binary = "std", thresholds = c(m = .1)))
  md <- match.data(mf); md$tr2 <- ifelse(md$treatment_group == cf[1], 1, 0)
  dd <- svydesign(ids = ~subclass, weights = ~weights, data = md)
  for (oc in c("VRR_binary","Comp_binary")) {
    co <- summary(svyglm(reformulate("tr2", oc), design = dd,
                         family = quasibinomial))$coefficients["tr2", ]
    cat(sprintf("  %-12s OR %.2f (%.2f-%.2f), p=%.4g\n", oc, exp(co[1]),
                exp(co[1]-1.96*co[2]), exp(co[1]+1.96*co[2]), co[4]))
  }
}

cat("\n--- S4. Sensitivity: strict threshold HR < 100 (exactly 100 counted as failure) ---\n")
df$VRR_strict <- ifelse(df$Pulse_Control < 100, 1, 0)
vs <- table(df$treatment_group, df$VRR_strict)[, c("1","0")]
print(vs)
cat(sprintf("Success: %s\n", paste(sprintf("%s %.1f%%", rownames(vs),
      100*vs[,1]/rowSums(vs)), collapse=" | ")))
for (cf in group_pairs) cat(or_wald(df, cf[1], cf[2], "VRR_strict"), "\n")

cat("\n--- S5. Cumulative dose and success (descriptive) ---\n")
print(aggregate(total_mgkg ~ treatment_group + VRR_binary, df,
                function(x) sprintf("%.3f (n=%d)", mean(x), length(x))))


## ---- 12. CLOSING ------------------------------------------------------
cat("\n", strrep("=",72), "\n### SESSION INFO\n", sep="")
print(sessionInfo())
cat("\nOUTPUT FOLDER:", OUT_DIR, "\n")
while (sink.number() > 0) sink()
cat("Done. Log:", file.path(OUT_DIR, "analysis_log_v2.txt"), "\n")
