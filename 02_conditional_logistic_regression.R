## =====================================================================
## DILTIAZEM FOR ATRIAL FIBRILLATION WITH RAPID VENTRICULAR RESPONSE
## Script 2 of 3  --  CONDITIONAL LOGISTIC REGRESSION ON THE MATCHED PAIRS
## =====================================================================
## Section 9a of script 01 reports an odds ratio from the matched 2x2
## table, which treats the matched observations as two independent
## groups. This script rebuilds exactly the same matching and then
## estimates the pair-conditional odds ratio, which respects the
## within-pair dependence that matching creates.
##
## It feeds the manuscript row in Table 5, panel A:
##   "Propensity score matching, 1:1, conditional logistic regression (ATT)"
##
## HOW TO RUN
##   Put this script in the same folder as "01_main_analysis.R" and run
##   it. If script 01 has already been run in this R session its objects
##   are reused; otherwise this script runs script 01 first, by itself.
##   Nothing has to be edited.
##
##   From a terminal:   Rscript 02_conditional_logistic_regression.R
##   From RStudio:      open the file and press "Source"
##
## OUTPUT
##   output_v2/clogit_analysis.log
##   output_v2/clogit_summary.csv
##
## Written for R >= 4.1. Tested with R 4.3.
## =====================================================================


## ---- 0. PORTABLE SETUP ----------------------------------------------
options(stringsAsFactors = FALSE, OutDec = ".", warn = 1)
invisible(suppressWarnings(Sys.setlocale("LC_COLLATE", "C")))
Sys.setenv(LANGUAGE = "en")

while (sink.number() > 0) sink()
options(error = function() {
  while (sink.number() > 0) sink()
  if (!interactive()) quit(save = "no", status = 1, runLast = FALSE)
})

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
HERE <- script_dir()

if (!nzchar(getOption("repos")[["CRAN"]]) ||
    identical(unname(getOption("repos")[["CRAN"]]), "@CRAN@"))
  options(repos = c(CRAN = "https://cloud.r-project.org"))

## No type = "binary" here: that argument is a hard error on Linux.
for (pk in c("MatchIt", "survival")) {
  if (!requireNamespace(pk, quietly = TRUE)) {
    message("Installing missing package: ", pk)
    try(utils::install.packages(pk), silent = TRUE)
  }
  if (!requireNamespace(pk, quietly = TRUE))
    stop("Package '", pk, "' is required but could not be installed.", call. = FALSE)
  suppressPackageStartupMessages(library(pk, character.only = TRUE))
}


## ---- 1. MAKE SURE THE MAIN ANALYSIS HAS BEEN RUN ---------------------
## exists()/get() search the whole search path by default. Base R already
## has stats::df (the F density), so inherits = FALSE is essential: without
## it, 'df' resolves to a FUNCTION whenever script 01 has not been run.
in_global  <- function(nm) exists(nm, envir = .GlobalEnv, inherits = FALSE)
get_global <- function(nm) get(nm, envir = .GlobalEnv, inherits = FALSE)

cohort_ready <- function() {
  if (!in_global("df")) return(FALSE)
  d <- get_global("df")
  is.data.frame(d) && nrow(d) == 214L &&
    all(c("treatment_group", "VRR_binary", "Comp_binary") %in% names(d)) &&
    nlevels(droplevels(factor(d$treatment_group))) == 3L
}

if (!cohort_ready()) {
  main_script <- file.path(HERE, "01_main_analysis.R")
  if (!file.exists(main_script))
    stop("Script 01 has not been run and '01_main_analysis.R' is not in\n",
         "  this folder (", HERE, ").\n",
         "  Run 01_main_analysis.R first, then run this script.", call. = FALSE)
  message("Script 01 has not been run in this session - running it now ...")
  source(main_script, local = FALSE, echo = FALSE)
  while (sink.number() > 0) sink()
  if (!cohort_ready())
    stop("Running 01_main_analysis.R did not produce the expected cohort.",
         call. = FALSE)
}


run_clogit <- function() {

  set.seed(2026)                       # same seed as script 01

  ## ---- 2. Output folder -----------------------------------------------
  out <- if (in_global("OUT_DIR")) get_global("OUT_DIR") else file.path(HERE, "output_v2")
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  log_path <- file.path(out, "clogit_analysis.log")

  ## ---- 3. Data ---------------------------------------------------------
  dd <- get_global("df")
  required <- c("treatment_group", "VRR_binary", "Comp_binary")
  if (length(setdiff(required, names(dd))))
    stop("Missing column(s): ", paste(setdiff(required, names(dd)), collapse = ", "),
         call. = FALSE)

  ## ---- 4. Reuse the objects built by script 01 -------------------------
  pairs_list <- if (in_global("group_pairs")) get_global("group_pairs") else
    list(c("Low_Dose","Standard_Bolus"), c("Low_Dose","Standard_Infusion"),
         c("Standard_Bolus","Standard_Infusion"))
  frm <- if (in_global("ps_f")) get_global("ps_f") else
    as.formula(paste("tr ~", paste(c("Age","Gender_bin","CAD_bin","New_Onset_AF_bin",
      "RateControl_bin","Pulse_Admission","Admission_SP","Admission_MAP","Weight"),
      collapse = " + ")))
  expected_pairs <- c(55, 44, 40)      # from analysis_log_v2.txt - a check only

  ## ---- 5. Log ----------------------------------------------------------
  while (sink.number() > 0) sink()
  sink(log_path, split = TRUE)
  on.exit(while (sink.number() > 0) sink(), add = TRUE)

  cat("CONDITIONAL LOGISTIC REGRESSION - CONDITIONAL ON THE MATCHED PAIRS\n")
  cat("Run on:", format(Sys.time()), "\n")
  cat("R:", R.version.string, "|", Sys.info()[["sysname"]], "\n")
  cat("Source: cohort from 01_main_analysis.R, n =", nrow(dd), "\n")
  cat("Matching: identical to section 9a (1:1 nearest neighbour,\n")
  cat("          caliper 0.2 SD on the logit scale, ATT)\n")
  cat("PS formula:", deparse(frm), "\n")
  cat("OR direction: the first group named is in the numerator (same as or_exact)\n")

  summary_tab <- NULL
  for (i in seq_along(pairs_list)) {
    cf <- pairs_list[[i]]
    cat("\n", paste(cf[1], "vs", cf[2]), "\n", sep = "")

    d2 <- droplevels(subset(dd, treatment_group %in% cf))
    d2$tr <- ifelse(d2$treatment_group == cf[2], 1, 0)   # same direction as 9a

    ps <- try(matchit(frm, data = d2, method = "nearest", distance = "glm",
                      link = "linear.logit", caliper = 0.2, replace = FALSE),
              silent = TRUE)
    if (inherits(ps, "try-error")) { cat("  Matching failed:\n"); print(ps); next }

    md <- match.data(ps)
    n_pairs <- sum(md$tr == 1)
    cat("  matched pairs:", n_pairs,
        if (n_pairs == expected_pairs[i]) "(same as analysis_log_v2)" else
          sprintf("  *** WARNING: analysis_log_v2 = %d ***", expected_pairs[i]), "\n")

    md$tr2      <- ifelse(md$treatment_group == cf[1], 1, 0)   # cf[1] in numerator
    md$subclass <- droplevels(factor(md$subclass))

    for (oc in c("VRR_binary", "Comp_binary")) {
      y  <- md[[oc]]
      tb <- table(group1 = md$tr2, outcome = y)
      cat("\n  [", oc, "]\n", sep = "")
      print(tb)

      ## A conditional likelihood takes information only from the pairs
      ## whose two members had different outcomes.
      discordant <- sum(tapply(y, md$subclass,
                              function(v) length(unique(v)) == 2L), na.rm = TRUE)
      cat("  discordant pairs:", discordant, "\n")

      if (!all(dim(tb) == c(2L, 2L)) || any(tb == 0) || discordant == 0) {
        cat("  -> zero cell / no discordant pairs: no conditional estimate possible\n")
        cat("     (the raw counts and the Firth estimate in section 9a remain valid)\n")
        next
      }

      fit <- try(clogit(as.formula(paste(oc, "~ tr2 + strata(subclass)")),
                        data = md), silent = TRUE)
      if (inherits(fit, "try-error")) { cat("  -> clogit returned an error\n"); next }

      s  <- summary(fit)
      or <- s$conf.int[1, 1]; lo <- s$conf.int[1, 3]; hi <- s$conf.int[1, 4]
      p  <- s$coefficients[1, 5]
      cat(sprintf("  -> OR %.2f (%.2f-%.2f), p=%.4g\n", or, lo, hi, p))
      summary_tab <- rbind(summary_tab, data.frame(
        comparison = paste(cf[1], "vs", cf[2]), outcome = oc,
        pairs = n_pairs, discordant = discordant, OR = round(or, 2),
        lo = round(lo, 2), hi = round(hi, 2), p = signif(p, 4),
        stringsAsFactors = FALSE))
    }
  }

  ## ---- 6. The Table 5 row ----------------------------------------------
  cat("\n", strrep("=", 72), "\n### SUMMARY\n", sep = "")
  if (!is.null(summary_tab)) {
    print(summary_tab, row.names = FALSE)
    write.csv(summary_tab, file.path(out, "clogit_summary.csv"), row.names = FALSE)
    v <- summary_tab[summary_tab$outcome == "VRR_binary", ]
    cat("\nTable 5 panel A -",
        "Propensity score matching, 1:1, conditional logistic regression (ATT):\n  ")
    cat(paste(sprintf("%.2f (%.2f-%.2f)", v$OR, v$lo, v$hi), collapse = " | "), "\n")
  } else cat("No estimate could be produced.\n")

  cat("\n", strrep("=", 72), "\n", sep = "")
  print(sessionInfo())
  cat("\nOUTPUT:", log_path, "\n")
  invisible(summary_tab)
}

clogit_summary <- run_clogit()
cat("Done. Log:", file.path(
  if (exists("OUT_DIR", envir = .GlobalEnv, inherits = FALSE))
    get("OUT_DIR", envir = .GlobalEnv, inherits = FALSE) else
    file.path(HERE, "output_v2"), "clogit_analysis.log"), "\n")
