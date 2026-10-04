# Diltiazem for atrial fibrillation with rapid ventricular response (RVR-AF)

Analysis code and de-identified dataset for the study comparing weight-based
low-dose diltiazem, standard-dose bolus and standard-dose infusion in the
emergency department.

Running the code reproduces every number in the manuscript, the supplement and
the tables, including the propensity-score analyses and the E-values.

The same code is provided in English (`english/`) and in Turkish (`turkish/`).
The two versions differ only in comments, messages and identifier names; they
read the same file and produce numerically identical results. Use whichever you
prefer — there is no need to run both.

---

## 1. Quick start

1. Download the whole deposit and keep the folder structure.
2. Install R (version 4.1 or later; developed and tested with R 4.3).
3. Run the scripts in order:

```
Rscript english/01_main_analysis.R
Rscript english/02_conditional_logistic_regression.R
Rscript english/03_evalue_verification.R
```

In RStudio, open each file and press **Source**. From an R console,
`source("english/01_main_analysis.R")`.

**Nothing has to be edited** — no file paths, no working directory, no encoding
setting. Each script locates its own folder, finds the data file, creates an
`output_v2` folder beside itself and writes everything there.

The scripts also work if you flatten the deposit and put all the files in one
folder, or if you move the data into a `data/` sub-folder.

### First run: packages

On the first run the scripts install any missing CRAN packages, so that run
needs an internet connection. Required:

| Package | Used for |
|---|---|
| `MatchIt` | propensity-score matching (nearest neighbour, optimal full) |
| `WeightIt` | stabilised IPTW and entropy balancing |
| `nnet` | the multinomial propensity model behind `WeightIt` |
| `cobalt` | balance diagnostics, love plots, overlap plots |
| `survey` | robust (design-based) standard errors for the weighted models |
| `tableone` | baseline table |
| `survival` | conditional logistic regression (script 02) |
| `EValue` | cross-check of the E-value calculation |
| `logistf` | Firth penalised logistic regression (sensitivity analyses) |

`EValue` and `logistf` are optional: if they cannot be installed, script 01
skips those two sensitivity analyses and says so in the log. Everything else is
required.

To install them yourself beforehand:

```r
install.packages(c("MatchIt", "WeightIt", "nnet", "cobalt", "survey",
                   "tableone", "survival", "EValue", "logistf"))
```

---

## 2. What is in the deposit

```
README_EN.md                     this file
README_TR.md                     the same information in Turkish
diltiazem_doz1set_eng_.csv       the de-identified dataset (n = 214)

english/
  01_main_analysis.R                     main analysis - run this first
  02_conditional_logistic_regression.R   pair-conditional estimate on the matched sample
  03_evalue_verification.R               independent re-derivation of the E-values

turkish/
  01_ana_analiz.R                        same as 01, in Turkish
  02_kosullu_lojistik_regresyon.R        same as 02, in Turkish
  03_evalue_dogrulama.R                  same as 03, in Turkish
```

### What each script does

**`01_main_analysis.R`** produces everything in one run: the dose descriptives
(Table 2), the baseline table and standardised mean differences (Table 1), the
primary outcome (Table 3), the secondary outcomes and adverse events (Table 4),
the infusion kinetics, the within-low-dose regression (Table 5), the
propensity-score matching and stabilised IPTW (Table 6), the E-values, and the
supplementary sensitivity analyses (BMI propensity model, entropy balancing,
optimal full matching, Firth correction, the strict `HR < 100` threshold).

Outputs: `output_v2/analysis_log_v2.txt` (full transcript),
`output_v2/table2_first_dose.csv`, `output_v2/table2b_cumulative_dose.csv`,
and the love/overlap plots as PNG.

**`02_conditional_logistic_regression.R`** rebuilds exactly the same 1:1
matching and fits a conditional logistic regression stratified on the matched
pair, which respects the within-pair dependence that matching creates. It feeds
the Table 5 panel A row *"Propensity score matching, 1:1, conditional logistic
regression (ATT)"*. If script 01 has not been run in the current R session,
this script runs it first by itself.

Outputs: `output_v2/clogit_analysis.log`, `output_v2/clogit_summary.csv`.

**`03_evalue_verification.R`** deliberately reuses nothing. It reads the raw
file again, rebuilds the nine covariates and the outcome from the raw columns,
and checks every derived quantity against the values printed in
`analysis_log_v2.txt` before reporting anything. Read the `VERDICT` line at the
end of its log before using any number from it.

Output: `output_v2/e_value.log` (plain ASCII).

---

## 3. Analysis summary

| | |
|---|---|
| Design | Prospective cohort, emergency department, n = 214 |
| Groups | Low dose 0.15 mg/kg (n = 105), standard bolus 0.25 mg/kg (n = 59), standard infusion 0.25 mg/kg over 10 min (n = 50) |
| Primary outcome | Heart rate ≤ 100 bpm, 10 minutes after the start of treatment |
| Primary adjusted analysis | Stabilised IPTW, ATE, body **weight** in the propensity model |
| Secondary adjusted analysis | 1:1 nearest-neighbour matching, caliper 0.2 SD on the logit scale, ATT |
| Propensity model covariates | age, sex, coronary artery disease, new-onset AF, chronic rate-control medication, admission heart rate, admission systolic pressure, admission mean arterial pressure, body weight |
| Unmeasured confounding | E-value from the weighted risk ratio (VanderWeele & Ding 2017) |

Two conventions are worth stating explicitly because they are easy to misread:

* In the matched analyses the **second** group named in a pair is treated as the
  "treated" group, so the estimand is the ATT with respect to that group.
* In every odds ratio and risk ratio the **first** group named is in the
  numerator.

The E-value is computed directly from an estimated risk ratio. The
`sqrt(OR)` approximation is **not** used, because the primary outcome is common
(55–90%) and that approximation is invalid in that range.

No p-value column is reported for the baseline table: balance is judged by
standardised mean differences.

---

## 4. Checking that your run reproduced the study

Script 01 stops with an explicit error if the dataset is not the archived one
(n = 214, group sizes 105 / 59 / 50). Beyond that, these are the values your
run should print:

| Quantity | Low dose | Standard bolus | Standard infusion |
|---|---|---|---|
| n | 105 | 59 | 50 |
| Primary outcome met | 58 (55.2%) | 53 (89.8%) | 44 (88.0%) |
| Mean age, years | 65.65 | 62.22 | 66.38 |
| Mean weight, kg | 74.57 | 75.53 | 75.42 |
| Mean admission heart rate | 149.22 | 148.61 | 153.04 |
| Mean admission systolic pressure | 130.19 | 131.36 | 131.20 |
| Mean admission MAP | 98.57 | 99.49 | 98.93 |

Omnibus chi-square for the primary outcome: X² = 30.56, df = 2, p = 2.32 × 10⁻⁷.
Matched pairs: 55, 44 and 40 for the three comparisons.
Script 03 prints `VERDICT: PASS` when its independent re-derivation agrees with
the published risk ratios.

Small differences in the last decimal place of the propensity-score results are
possible across package versions; the group sizes, event counts and descriptive
statistics above are exact and should match byte for byte.

---

## 5. Dataset

`diltiazem_doz1set_eng_.csv` — 214 rows, 52 columns, one row per patient.
Semicolon-separated, **comma as the decimal mark**, CRLF line endings. The file
is ASCII apart from a single byte: a Turkish letter inside the column name
`@yari_doz`, which no analysis uses.

The columns the analysis uses:

| Column | Meaning | Coding |
|---|---|---|
| `Patient_NB` | patient number | 1–214 |
| `treatment_group` | treatment arm | 1 = standard infusion, 2 = standard bolus, 3 = low dose |
| `Age` | age | years |
| `Gender` | sex | 1 = male |
| `Weight`, `Height`, `BMI` | body size | kg, cm, kg/m² |
| `Ordered_Dose` | first diltiazem dose actually given | mg (pump-delivered, no rounding) |
| `CAD` | coronary artery disease | 1 = yes |
| `New_Onset_AF` | new-onset atrial fibrillation | 1 = yes |
| `Use_of_Antiarrhytmic` | chronic rate-control medication | 1 = yes |
| `Which_antiarrhytmic` | class of that medication | 1 = calcium-channel blocker, 2 = beta blocker, 3 = none, 4 = digoxin |
| `Pulse_Admission` | heart rate on arrival | bpm |
| `Pulse_Control` | heart rate at 10 minutes | bpm |
| `First_Dose_Ventricular_Rate_Response` | recorded primary outcome | 1 = rate controlled, 2 = not controlled |
| `Admission_SP`, `Admission_DP`, `Admission_MAP` | blood pressure on arrival | mmHg |
| `Control_Systolic_Blood_Pressure`, `Control_Diastole`, `Control_MAP` | blood pressure at 10 minutes | mmHg |
| `Complication` | adverse event | 1 = yes |
| `Result` | disposition | 2 = admitted to hospital |
| `@150sec_Pulse` … `@600sec_Pulse` | heart rate at 150/300/450/600 s | bpm, infusion group only |
| `@150sec_VRR` … `@600sec_VRR` | rate controlled at that time point | 1 = yes, infusion group only |
| `Low_dose`, `aditional_dose`, `additional_medication`, `@0.10_mg_dose`, `@0.35_mg_dose`, `@Beta_Bloker_Ventricular_Rate_Response` | rescue treatment given after the first dose | see the cumulative-exposure block in script 01 |
| `Secound_Dose_Ventricular_Rate_Response` | response to the rescue dose | |
| `After_First_Dose_HR`, `After_Sec_Dose_Hr`, `After_BB_HR`, `LD_0,10_HR`, `LD_0,35_HR` | heart rate after the respective dose | bpm |

The code derives the binary analysis variables from these columns; it never
reads a pre-computed indicator. The data file is left exactly as archived.

Two integrity checks run automatically before any analysis: the recorded
primary-outcome flag is compared patient by patient against the
`heart rate ≤ 100 bpm` rule, and the ordered dose is compared against
`body weight × protocol coefficient`. Both must agree for every patient.

---

## 6. Portability notes

These notes explain why the code is written the way it is. They matter if you
modify it.

**The data file is read as raw bytes, not through an encoding connection.**
The original working scripts passed `fileEncoding = "latin1"` to `read.csv()`.
That works on a Turkish Windows machine, but on Linux R aborts the re-encoding
connection at the single non-ASCII byte in the header and returns 19 columns
and 0 rows (verified with R 4.3.3). The scripts therefore read the file as
bytes, replace any non-ASCII byte with `.` — which is what `make.names()` would
do to it anyway — and parse the text from memory. Identical bytes then give an
identical data frame on every operating system and in every locale.

**Collation is fixed to the C locale.** `Sys.setlocale("LC_COLLATE", "C")` is
set at the top of every script, so `sort()` and factor level order do not depend
on the system language.

**`install.packages()` is called without `type = "binary"`.** That argument is a
hard error on Linux and is unnecessary on Windows and macOS.

**Paths are resolved at run time.** Each script determines its own folder
whether it was started with `Rscript`, with `source()`, or with the RStudio
*Source* button, and then searches for the data file beside itself, one level
up, in a `data/` sub-folder and in the working directory. No `setwd()` to a
fixed path.

**Plots degrade gracefully.** Figures are written with `ggplot2::ggsave()` and
fall back to `png()`; on a headless server with no graphics device the script
notes that the figure was skipped and carries on.

**The script files are pure ASCII**, including the Turkish versions, which use
diacritic-free Turkish. This is deliberate: the comments then read correctly
whatever locale the file is opened in.

**If a run is interrupted** while the log file is open, the next run closes the
stray connection by itself. If your R console still seems to swallow output
after a crash, run `while (sink.number() > 0) sink()` once.

---

## 7. Reuse

The dataset is de-identified: it contains no names, no dates, no identifiers and
no free text. Ethics approval and consent are described in the manuscript.

If you reuse this code or dataset, please cite the published article and this
deposit.

<!-- Add the article citation, the Zenodo DOI and your chosen licence here
     before publishing the record. -->
