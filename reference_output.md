# Reference output

Put the logs from the authors' own run here, so that anyone who re-runs the
code can compare their results line by line against the published ones:

    analysis_log_v2.txt     from english/01_main_analysis.R
    clogit_analysis.log     from english/02_conditional_logistic_regression.R
    e_value.log             from english/03_evalue_verification.R

A fresh run writes its own copies into `english/output_v2/` (or
`turkish/cikti_v2/`), which git ignores. Comparing the two is the quickest
reproducibility check there is: every number should be identical.

Before copying the logs here, re-run the scripts from a neutral folder
(for example `C:\diltiazem_zenodo\`). The logs record the full path of the
folder they ran in, and that path should not contain a personal user name or
a private folder name once the record is public.

### DOSE ANALYSES

-- Ordered_Dose check (body weight x protocol coefficient) --
Exact agreement: 214/214 | largest deviation 0.0000 mg

-- TABLE 2: first dose --
                    n               weight                      dose_mg           mg_kg
Low_Dose          105 74.6 +/- 7.8 (57-96)  11.19 +/- 1.18 (8.55-14.40) 0.15 (protocol)
Standard_Bolus     59 75.5 +/- 8.3 (58-96) 18.88 +/- 2.08 (14.50-24.00) 0.25 (protocol)
Standard_Infusion  50 75.4 +/- 7.5 (58-88) 18.86 +/- 1.87 (14.50-22.00) 0.25 (protocol)

-- Rescue treatment --
                   extra dose (1=yes)
group                1  2
  Low_Dose          43 62
  Standard_Bolus     6 53
  Standard_Infusion  5 45

0.35 mg/kg rescue:
                   0.35
group                0  1  2
  Low_Dose          88 11  6
  Standard_Bolus    53  3  3
  Standard_Infusion 45  4  1

Beta blocker:
                   BB
group                0  1  2
  Low_Dose          99  3  3
  Standard_Bolus    56  2  1
  Standard_Infusion 49  0  1

-- CUMULATIVE EXPOSURE (first dose + rescue) --
                    n first_mgkg      total_mgkg      total_mg received_extra
Low_Dose          105       0.15 0.248 +/- 0.161 18.5 +/- 12.2             43
Standard_Bolus     59       0.25 0.286 +/- 0.107  21.6 +/- 8.4              6
Standard_Infusion  50       0.25 0.285 +/- 0.106  21.7 +/- 9.3              5

Cumulative mg/kg across groups (Welch):
  F=2.013 df=(2, 128.8) p=0.1378

========================================================================
### TABLE 1 - BASELINE
                             Stratified by treatment_group
                              Low_Dose       Standard_Bolus Standard_Infusion
  n                              105             59             50           
  Age (mean (SD))              65.65 (11.50)  62.22 (16.03)  66.38 (13.56)   
  Gender_bin = 1 (%)              35 (33.3)      26 (44.1)      23 (46.0)    
  CAD_bin = 1 (%)                  9 ( 8.6)       5 ( 8.5)       7 (14.0)    
  New_Onset_AF_bin = 1 (%)        27 (25.7)      12 (20.3)      14 (28.0)    
  RateControl_bin = 1 (%)         80 (76.2)      48 (81.4)      37 (74.0)    
  Weight (mean (SD))           74.57 (7.84)   75.53 (8.32)   75.42 (7.48)    
  BMI (mean (SD))              25.62 (1.96)   26.24 (2.29)   26.36 (2.01)    
  Pulse_Admission (mean (SD)) 149.22 (11.19) 148.61 (11.66) 153.04 (10.99)   
  Admission_SP (mean (SD))    130.19 (14.93) 131.36 (13.83) 131.20 (14.38)   
  Admission_MAP (mean (SD))    98.57 (8.48)   99.49 (7.45)   98.93 (8.82)    

-- SMD before adjustment (largest of the pairwise comparisons, binary='std') --
Balance summary across all treatment pairs
                    Type Max.Diff.Un     M.Threshold.Un
Age              Contin.      0.3010 Not Balanced, >0.1
Gender_bin        Binary      0.2591 Not Balanced, >0.1
CAD_bin           Binary      0.1821 Not Balanced, >0.1
New_Onset_AF_bin  Binary      0.1782 Not Balanced, >0.1
RateControl_bin   Binary      0.1758 Not Balanced, >0.1
Weight           Contin.      0.1210 Not Balanced, >0.1
BMI              Contin.      0.3518 Not Balanced, >0.1
Pulse_Admission  Contin.      0.3926 Not Balanced, >0.1
Admission_SP     Contin.      0.0810     Balanced, <0.1
Admission_MAP    Contin.      0.1112 Not Balanced, >0.1

Balance tally for mean differences
                   count
Balanced, <0.1         1
Not Balanced, >0.1     9

Variable with the greatest mean difference
        Variable Max.Diff.Un     M.Threshold.Un
 Pulse_Admission      0.3926 Not Balanced, >0.1

Sample sizes
    Low_Dose Standard_Bolus Standard_Infusion
All      105             59                50

-- Class of chronic rate-control drug --
                   
                    CCB Beta_blocker Digoxin None
  Low_Dose           61           19       0   25
  Standard_Bolus     33           14       1   11
  Standard_Infusion  23           14       0   13

BMI range: 20.2 - 31.6 kg/m2 | BMI >= 30: 4 patients

========================================================================
### TABLE 3 - PRIMARY OUTCOME (HR <= 100 bpm at 10 min)
                   
                     1  0
  Low_Dose          58 47
  Standard_Bolus    53  6
  Standard_Infusion 44  6
Success: Low_Dose 55.2% | Standard_Bolus 89.8% | Standard_Infusion 88.0%
Omnibus chi-square: X2=30.56, df=2, p=2.32e-07

Pairwise comparisons, Bonferroni corrected (chi-square with Yates correction):

	Pairwise comparisons using Pairwise comparison of proportions 

data:  vrr[, 1] out of rowSums(vrr) 

                  Low_Dose Standard_Bolus
Standard_Bolus    3.7e-05  -             
Standard_Infusion 0.00037  1.00000       

P value adjustment method: bonferroni 

-- Unadjusted OR (Wald) --
Low_Dose vs Standard_Bolus             OR 0.14 (0.06-0.35), p=3.205e-05 
Low_Dose vs Standard_Infusion          OR 0.17 (0.07-0.43), p=0.0001892 
Standard_Bolus vs Standard_Infusion    OR 1.20 (0.36-4.00), p=0.7612 

-- Crude risk difference --
Low_Dose vs Standard_Bolus             RD -34.6% (-46.8 to -22.3) 
Low_Dose vs Standard_Infusion          RD -32.8% (-45.9 to -19.7) 
Standard_Bolus vs Standard_Infusion    RD +1.8% (-10.0 to 13.7) 

========================================================================
### TABLE 4 - SECONDARY OUTCOMES
Pulse_Control                     108.9 +/- 27.2 | 85.1 +/- 19.5 | 87.0 +/- 23.2
                                  Welch F=24.095 df=(2,121.2) p=1.546e-09 | Bartlett p=0.02028
Control_Systolic_Blood_Pressure   121.6 +/- 13.3 | 115.6 +/- 16.2 | 119.0 +/- 17.5
                                  Welch F=2.998 df=(2,102.8) p=0.05426 | Bartlett p=0.04936
Control_MAP                       92.5 +/- 9.0 | 88.4 +/- 13.8 | 90.1 +/- 14.1
                                  Welch F=2.374 df=(2,94.7) p=0.09863 | Bartlett p=5.034e-05

-- Adverse events --
                   
                      1   0
  Low_Dose            2 103
  Standard_Bolus      5  54
  Standard_Infusion   4  46
Omnibus Fisher p = 0.0924 
Low_Dose vs Standard_Bolus             OR 0.21 (0.04-1.12), p=0.06717 
Low_Dose vs Standard_Infusion          OR 0.22 (0.04-1.26), p=0.08987 
Standard_Bolus vs Standard_Infusion    OR 1.06 (0.27-4.20), p=0.9285 

-- Breakdown by event type --
 Patient_NB   treatment_group Pulse_Control Control_Systolic_Blood_Pressure hypotension bradycardia
          9 Standard_Infusion            84                              60        TRUE       FALSE
         13 Standard_Infusion           100                              80        TRUE       FALSE
         24 Standard_Infusion           130                              70        TRUE       FALSE
         30 Standard_Infusion            80                              80        TRUE       FALSE
         61    Standard_Bolus            76                              60        TRUE       FALSE
         68    Standard_Bolus            90                              70        TRUE       FALSE
         72    Standard_Bolus            76                              80        TRUE       FALSE
         78    Standard_Bolus            74                              70        TRUE       FALSE
         79    Standard_Bolus            82                              90       FALSE       FALSE
        121          Low_Dose           130                              80        TRUE       FALSE
        156          Low_Dose            88                              70        TRUE       FALSE
Total 11 events | hypotension 10 | bradycardia 0 | other 1

-- Hospital admission (post hoc) --
                   
                     1  0
  Low_Dose          11 94
  Standard_Bolus     7 52
  Standard_Infusion  8 42
Fisher p = 0.5862 

========================================================================
### INFUSION KINETICS (standard infusion, n=50)
Source: the recorded response columns (not recomputed from heart rate)

 150 s: response column  4/50 =  8% | from heart rate HR<=100:  4 
 300 s: response column 23/50 = 46% | from heart rate HR<=100: 23 
 450 s: response column 39/50 = 78% | from heart rate HR<=100: 39 
 600 s: response column 44/50 = 88% | from heart rate HR<=100: 44 

Paired t-test against baseline, Bonferroni (4 comparisons):
   150 s: p=4.47e-05 (Bonferroni 0.000179)
   300 s: p=1.25e-15 (Bonferroni 4.98e-15)
   450 s: p=7.32e-22 (Bonferroni 2.93e-21)
   600 s: p=5.21e-24 (Bonferroni 2.09e-23)

========================================================================
### TABLE 5 - PREDICTORS WITHIN THE LOW-DOSE GROUP
                    OR    lo    hi     p
Age              1.030 0.990 1.073 0.143
Gender_bin       1.416 0.607 3.304 0.421
CAD_bin          0.767 0.177 3.321 0.723
New_Onset_AF_bin 0.994 0.177 5.582 0.995
RateControl_bin  0.873 0.151 5.053 0.879
n=105, successes=58, failures=47 (EPV = 9.4)
Model fit: LR chi-square = 2.96, df = 5, p = 0.707
Firth sensitivity analysis:
                    OR    lo    hi
Age              1.028 0.990 1.070
Gender_bin       1.383 0.613 3.188
CAD_bin          0.767 0.190 3.147
New_Onset_AF_bin 0.986 0.192 5.140
RateControl_bin  0.884 0.164 4.617

========================================================================
### TABLE 6 - PROPENSITY SCORE ANALYSES

--- 9a. Matching (1:1 nearest neighbour, caliper 0.2 SD on the logit scale, ATT) ---
    In each pair the SECOND group named is treated as the 'treated'
    group, so the estimand is the ATT with respect to that group.
    Odds ratios are reported with the FIRST group in the numerator.

 Low_Dose vs Standard_Bolus - matched pairs: 55 
Balance Measures
                     Type Diff.Un Diff.Adj        M.Threshold
distance         Distance  0.3978   0.0156     Balanced, <0.1
Age               Contin. -0.2138   0.0544     Balanced, <0.1
Gender_bin         Binary  0.2162  -0.1831 Not Balanced, >0.1
CAD_bin            Binary -0.0035   0.0653     Balanced, <0.1
New_Onset_AF_bin   Binary -0.1335  -0.1355 Not Balanced, >0.1
RateControl_bin    Binary  0.1326   0.1867 Not Balanced, >0.1
Pulse_Admission   Contin. -0.0522  -0.1716 Not Balanced, >0.1
Admission_SP      Contin.  0.0843  -0.0657     Balanced, <0.1
Admission_MAP     Contin.  0.1235  -0.0407     Balanced, <0.1
Weight            Contin.  0.1146  -0.0459     Balanced, <0.1

Balance tally for mean differences
                   count
Balanced, <0.1         6
Not Balanced, >0.1     4

Variable with the greatest mean difference
        Variable Diff.Adj        M.Threshold
 RateControl_bin   0.1867 Not Balanced, >0.1

Sample sizes
          Control Treated
All           105      59
Matched        55      55
Unmatched      50       4
  >> Covariates above the 0.1 threshold: 4 / 9
  VRR_binary   OR 0.21 (0.06-0.65), p=0.004196
  Comp_binary  OR 0.24 (0.00-2.52), p=0.3634

 Low_Dose vs Standard_Infusion - matched pairs: 44 
Balance Measures
                     Type Diff.Un Diff.Adj        M.Threshold
distance         Distance  0.5345   0.0206     Balanced, <0.1
Age               Contin.  0.0540   0.0251     Balanced, <0.1
Gender_bin         Binary  0.2541   0.0000     Balanced, <0.1
CAD_bin            Binary  0.1564   0.0000     Balanced, <0.1
New_Onset_AF_bin   Binary  0.0509   0.2025 Not Balanced, >0.1
RateControl_bin    Binary -0.0499  -0.1554 Not Balanced, >0.1
Pulse_Admission   Contin.  0.3476   0.0000     Balanced, <0.1
Admission_SP      Contin.  0.0702   0.0632     Balanced, <0.1
Admission_MAP     Contin.  0.0410   0.0859     Balanced, <0.1
Weight            Contin.  0.1135  -0.2037 Not Balanced, >0.1

Balance tally for mean differences
                   count
Balanced, <0.1         7
Not Balanced, >0.1     3

Variable with the greatest mean difference
 Variable Diff.Adj        M.Threshold
   Weight  -0.2037 Not Balanced, >0.1

Sample sizes
          Control Treated
All           105      50
Matched        44      44
Unmatched      61       6
  >> Covariates above the 0.1 threshold: 3 / 9
  VRR_binary   OR 0.13 (0.04-0.41), p=9.486e-05
  Comp_binary  0/44 vs 4/44 - NOT ESTIMABLE (zero cell)
               Firth OR 0.10 (0.00-0.99)

 Standard_Bolus vs Standard_Infusion - matched pairs: 40 
Balance Measures
                     Type Diff.Un Diff.Adj        M.Threshold
distance         Distance  0.7247   0.0486     Balanced, <0.1
Age               Contin.  0.3068   0.0793     Balanced, <0.1
Gender_bin         Binary  0.0388  -0.0502     Balanced, <0.1
CAD_bin            Binary  0.1592   0.0000     Balanced, <0.1
New_Onset_AF_bin   Binary  0.1706   0.1114 Not Balanced, >0.1
RateControl_bin    Binary -0.1677  -0.0570     Balanced, <0.1
Pulse_Admission   Contin.  0.4030  -0.0637     Balanced, <0.1
Admission_SP      Contin. -0.0108   0.0000     Balanced, <0.1
Admission_MAP     Contin. -0.0633   0.0000     Balanced, <0.1
Weight            Contin. -0.0141  -0.1338 Not Balanced, >0.1

Balance tally for mean differences
                   count
Balanced, <0.1         8
Not Balanced, >0.1     2

Variable with the greatest mean difference
 Variable Diff.Adj        M.Threshold
   Weight  -0.1338 Not Balanced, >0.1

Sample sizes
          Control Treated
All            59      50
Matched        40      40
Unmatched      19      10
  >> Covariates above the 0.1 threshold: 2 / 9
  VRR_binary   OR 1.00 (0.21-4.77), p=1
  Comp_binary  OR 0.48 (0.04-3.57), p=0.6752

--- 9b. IPTW (stabilised, ATE) - PRIMARY ADJUSTED ANALYSIS ---
                  Summary of weights

- Weight ranges:

                    Min                                 Max
Low_Dose          0.76       |----------|             1.546
Standard_Bolus    0.396 |---------------------------| 2.385
Standard_Infusion 0.48   |----------------------|     2.106

- Units with the 5 most extreme weights by group:
                                                
                     214   180   179   155   153
          Low_Dose 1.378 1.406 1.424 1.475 1.546
                      97    52    71    92    57
    Standard_Bolus 1.731 1.889   2.1 2.213 2.385
                      25    45    37    14    24
 Standard_Infusion 1.655 1.922 1.986 2.081 2.106

- Weight statistics:

                  Coef of Var   MAD Entropy # Zeros
Low_Dose                0.172 0.138   0.014       0
Standard_Bolus          0.424 0.326   0.081       0
Standard_Infusion       0.432 0.34    0.084       0

- Effective Sample Sizes:

           Low_Dose Standard_Bolus Standard_Infusion
Unweighted   105.            59.               50.  
Weighted     102.02          50.14             42.26
Balance summary across all treatment pairs
                    Type Max.Diff.Un Max.Diff.Adj    M.Threshold
Age              Contin.      0.3010       0.0851 Balanced, <0.1
Gender_bin        Binary      0.2591       0.0355 Balanced, <0.1
CAD_bin           Binary      0.1821       0.0291 Balanced, <0.1
New_Onset_AF_bin  Binary      0.1782       0.0210 Balanced, <0.1
RateControl_bin   Binary      0.1758       0.0436 Balanced, <0.1
Pulse_Admission  Contin.      0.3926       0.0263 Balanced, <0.1
Admission_SP     Contin.      0.0810       0.0543 Balanced, <0.1
Admission_MAP    Contin.      0.1112       0.0637 Balanced, <0.1
Weight           Contin.      0.1210       0.0362 Balanced, <0.1

Balance tally for mean differences
                   count
Balanced, <0.1         9
Not Balanced, >0.1     0

Variable with the greatest mean difference
 Variable Max.Diff.Adj    M.Threshold
      Age       0.0851 Balanced, <0.1

Effective sample sizes
           Low_Dose Standard_Bolus Standard_Infusion
Unadjusted   105.            59.               50.  
Adjusted     102.02          50.14             42.26

IPTW odds ratios (robust SE):
 -- VRR_binary --
   Low_Dose vs Standard_Bolus             OR 0.13 (0.05-0.38), p=0.0001771 
   Low_Dose vs Standard_Infusion          OR 0.17 (0.06-0.47), p=0.0009705 
   Standard_Bolus vs Standard_Infusion    OR 1.23 (0.32-4.77), p=0.7678 
 -- Comp_binary --
   Low_Dose vs Standard_Bolus             OR 0.32 (0.06-1.72), p=0.1849 
   Low_Dose vs Standard_Infusion          OR 0.15 (0.03-0.91), p=0.04021 
   Standard_Bolus vs Standard_Infusion    OR 0.48 (0.11-2.04), p=0.3206 

IPTW risk ratios (log-Poisson, robust SE) - primary outcome:
  Low_Dose vs Standard_Bolus             RR 0.62 (0.51-0.75), p=3.043e-06
  Low_Dose vs Standard_Infusion          RR 0.63 (0.51-0.78), p=2.238e-05
  Standard_Bolus vs Standard_Infusion    RR 1.02 (0.88-1.18), p=0.7699

========================================================================
### E-VALUE (from the weighted risk ratio)
Method: VanderWeele & Ding 2017. The risk ratio was estimated directly;
the sqrt(OR) approximation was NOT used (it is invalid for a common outcome).

Low_Dose vs Standard_Bolus:
  RR 0.617 (0.507-0.750) -> E-value 2.62 (at the CI limit 2.00)
  Cross-check with the EValue package:
             point     lower     upper
RR       0.6168009 0.5071233 0.7501989
E-values 2.6248839        NA 1.9992045
Low_Dose vs Standard_Infusion:
  RR 0.630 (0.513-0.775) -> E-value 2.55 (at the CI limit 1.90)
  Cross-check with the EValue package:
             point     lower     upper
RR       0.6303746 0.5126439 0.7751427
E-values 2.5508138        NA 1.9018321

========================================================================
### SUPPLEMENT - SENSITIVITY ANALYSES

--- S1. IPTW with BMI in the PS model instead of weight ---
Balance summary across all treatment pairs
                    Type Max.Diff.Adj        M.Threshold
Age              Contin.       0.0696     Balanced, <0.1
Gender_bin        Binary       0.0470     Balanced, <0.1
CAD_bin           Binary       0.0201     Balanced, <0.1
New_Onset_AF_bin  Binary       0.0356     Balanced, <0.1
RateControl_bin   Binary       0.0510     Balanced, <0.1
Pulse_Admission  Contin.       0.0320     Balanced, <0.1
Admission_SP     Contin.       0.0675     Balanced, <0.1
Admission_MAP    Contin.       0.0753     Balanced, <0.1
BMI              Contin.       0.1211 Not Balanced, >0.1

Balance tally for mean differences
                   count
Balanced, <0.1         8
Not Balanced, >0.1     1

Variable with the greatest mean difference
 Variable Max.Diff.Adj        M.Threshold
      BMI       0.1211 Not Balanced, >0.1

Effective sample sizes
           Low_Dose Standard_Bolus Standard_Infusion
Unadjusted   105.             59.              50.  
Adjusted      98.99           49.5             42.52
 -- VRR_binary --
   Low_Dose vs Standard_Bolus             OR 0.13 (0.05-0.38), p=0.000222 
   Low_Dose vs Standard_Infusion          OR 0.13 (0.04-0.36), p=0.0001469 
   Standard_Bolus vs Standard_Infusion    OR 0.95 (0.24-3.75), p=0.9448 
 -- Comp_binary --
   Low_Dose vs Standard_Bolus             OR 0.25 (0.05-1.39), p=0.116 
   Low_Dose vs Standard_Infusion          OR 0.14 (0.02-0.81), p=0.03002 
   Standard_Bolus vs Standard_Infusion    OR 0.54 (0.13-2.30), p=0.4078 

--- S2. Entropy balancing (weight) ---
                  Summary of weights

- Weight ranges:

                    Min                                  Max
Low_Dose          1.437 |----|                         3.272
Standard_Bolus    1.444 |-----------------|            7.155
Standard_Infusion 1.408 |---------------------------| 10.158

- Units with the 5 most extreme weights by group:
                                                 
                     192   166   150   180    179
          Low_Dose 2.665 2.856 2.893 2.984  3.272
                      57    97    52    92     71
    Standard_Bolus 5.359 5.636 6.169 6.917  7.155
                      25    14    37    45     24
 Standard_Infusion 7.679 7.689 8.744 8.746 10.158

- Weight statistics:

                  Coef of Var   MAD Entropy # Zeros
Low_Dose                0.18  0.145   0.015       0
Standard_Bolus          0.345 0.272   0.057       0
Standard_Infusion       0.47  0.371   0.101       0

- Effective Sample Sizes:

           Low_Dose Standard_Bolus Standard_Infusion
Unweighted   105.            59.               50.  
Weighted     101.74          52.83             41.09
Balance summary across all treatment pairs
                    Type Max.Diff.Adj    M.Threshold
Age              Contin.            0 Balanced, <0.1
Gender_bin        Binary            0 Balanced, <0.1
CAD_bin           Binary            0 Balanced, <0.1
New_Onset_AF_bin  Binary            0 Balanced, <0.1
RateControl_bin   Binary            0 Balanced, <0.1
Pulse_Admission  Contin.            0 Balanced, <0.1
Admission_SP     Contin.            0 Balanced, <0.1
Admission_MAP    Contin.            0 Balanced, <0.1
Weight           Contin.            0 Balanced, <0.1

Balance tally for mean differences
                   count
Balanced, <0.1         9
Not Balanced, >0.1     0

Variable with the greatest mean difference
 Variable Max.Diff.Adj    M.Threshold
  CAD_bin            0 Balanced, <0.1

Effective sample sizes
           Low_Dose Standard_Bolus Standard_Infusion
Unadjusted   105.            59.               50.  
Adjusted     101.74          52.83             41.09
 -- VRR_binary --
   Low_Dose vs Standard_Bolus             OR 0.14 (0.05-0.37), p=0.0001407 
   Low_Dose vs Standard_Infusion          OR 0.16 (0.05-0.49), p=0.001447 
   Standard_Bolus vs Standard_Infusion    OR 1.20 (0.30-4.74), p=0.797 
 -- Comp_binary --
   Low_Dose vs Standard_Bolus             OR 0.30 (0.05-1.62), p=0.1624 
   Low_Dose vs Standard_Infusion          OR 0.15 (0.02-0.91), p=0.04037 
   Standard_Bolus vs Standard_Infusion    OR 0.50 (0.11-2.20), p=0.3607 

--- S3. Optimal full matching ---

 Low_Dose vs Standard_Bolus 
Balance Measures
                     Type Diff.Adj        M.Threshold
distance         Distance   0.0617     Balanced, <0.1
Age               Contin.  -0.2185 Not Balanced, >0.1
Gender_bin         Binary  -0.0589     Balanced, <0.1
CAD_bin            Binary  -0.0071     Balanced, <0.1
New_Onset_AF_bin   Binary  -0.0404     Balanced, <0.1
RateControl_bin    Binary   0.0823     Balanced, <0.1
Pulse_Admission   Contin.   0.0570     Balanced, <0.1
Admission_SP      Contin.  -0.2241 Not Balanced, >0.1
Admission_MAP     Contin.  -0.3497 Not Balanced, >0.1
Weight            Contin.  -0.0613     Balanced, <0.1

Balance tally for mean differences
                   count
Balanced, <0.1         7
Not Balanced, >0.1     3

Variable with the greatest mean difference
      Variable Diff.Adj        M.Threshold
 Admission_MAP  -0.3497 Not Balanced, >0.1

Sample sizes
                     Control Treated
All                   105.        59
Matched (ESS)          41.83      59
Matched (Unweighted)  105.        59
  VRR_binary   OR 0.18 (0.04-0.80), p=0.02949
  Comp_binary  OR 0.21 (0.03-1.63), p=0.1418

 Low_Dose vs Standard_Infusion 
Balance Measures
                     Type Diff.Adj        M.Threshold
distance         Distance  -0.0145     Balanced, <0.1
Age               Contin.   0.0514     Balanced, <0.1
Gender_bin         Binary   0.1712 Not Balanced, >0.1
CAD_bin            Binary   0.2133 Not Balanced, >0.1
New_Onset_AF_bin   Binary   0.0767     Balanced, <0.1
RateControl_bin    Binary  -0.0151     Balanced, <0.1
Pulse_Admission   Contin.  -0.2711 Not Balanced, >0.1
Admission_SP      Contin.   0.0709     Balanced, <0.1
Admission_MAP     Contin.   0.0758     Balanced, <0.1
Weight            Contin.   0.0882     Balanced, <0.1

Balance tally for mean differences
                   count
Balanced, <0.1         7
Not Balanced, >0.1     3

Variable with the greatest mean difference
        Variable Diff.Adj        M.Threshold
 Pulse_Admission  -0.2711 Not Balanced, >0.1

Sample sizes
                     Control Treated
All                   105.        50
Matched (ESS)          33.94      50
Matched (Unweighted)  105.        50
  VRR_binary   OR 0.12 (0.04-0.37), p=0.0007561
  Comp_binary  OR 0.15 (0.02-1.05), p=0.06318

 Standard_Bolus vs Standard_Infusion 
Balance Measures
                     Type Diff.Adj        M.Threshold
distance         Distance   0.0190     Balanced, <0.1
Age               Contin.  -0.1056 Not Balanced, >0.1
Gender_bin         Binary   0.0040     Balanced, <0.1
CAD_bin            Binary  -0.1489 Not Balanced, >0.1
New_Onset_AF_bin   Binary   0.1700 Not Balanced, >0.1
RateControl_bin    Binary  -0.1436 Not Balanced, >0.1
Pulse_Admission   Contin.   0.0637     Balanced, <0.1
Admission_SP      Contin.  -0.0359     Balanced, <0.1
Admission_MAP     Contin.  -0.0352     Balanced, <0.1
Weight            Contin.  -0.1860 Not Balanced, >0.1

Balance tally for mean differences
                   count
Balanced, <0.1         5
Not Balanced, >0.1     5

Variable with the greatest mean difference
 Variable Diff.Adj        M.Threshold
   Weight   -0.186 Not Balanced, >0.1

Sample sizes
                     Control Treated
All                    59.        50
Matched (ESS)          24.41      50
Matched (Unweighted)   59.        50
  VRR_binary   OR 1.79 (0.44-7.37), p=0.4235
  Comp_binary  OR 0.55 (0.10-3.03), p=0.5003

--- S4. Sensitivity: strict threshold HR < 100 (exactly 100 counted as failure) ---
                   
                     1  0
  Low_Dose          58 47
  Standard_Bolus    53  6
  Standard_Infusion 43  7
Success: Low_Dose 55.2% | Standard_Bolus 89.8% | Standard_Infusion 86.0%
Low_Dose vs Standard_Bolus             OR 0.14 (0.06-0.35), p=3.205e-05 
Low_Dose vs Standard_Infusion          OR 0.20 (0.08-0.49), p=0.0003881 
Standard_Bolus vs Standard_Infusion    OR 1.44 (0.45-4.60), p=0.5402 

--- S5. Cumulative dose and success (descriptive) ---
    treatment_group VRR_binary   total_mgkg
1          Low_Dose          0 0.368 (n=47)
2    Standard_Bolus          0  0.600 (n=6)
3 Standard_Infusion          0  0.542 (n=6)
4          Low_Dose          1 0.150 (n=58)
5    Standard_Bolus          1 0.250 (n=53)
6 Standard_Infusion          1 0.250 (n=44)

========================================================================
### SESSION INFO
R version 4.6.0 (2026-04-24 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] C
system code page: 65001

time zone: Europe/Istanbul
tzcode source: internal

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] EValue_4.1.4    logistf_1.26.1  nnet_7.3-20     tableone_0.13.2 survey_4.5      survival_3.8-6  Matrix_1.7-5   
 [8] cobalt_4.6.2    WeightIt_1.7.0  MatchIt_4.7.2  

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.1       dplyr_1.2.1            farver_2.1.2           S7_0.2.2               mathjaxr_2.0-0        
 [6] labelled_2.16.0        digest_0.6.39          rpart_4.1.27           lifecycle_1.0.5        magrittr_2.0.5        
[11] compiler_4.6.0         rlang_1.2.0            tools_4.6.0            labeling_0.4.3         RColorBrewer_1.1-3    
[16] withr_3.0.2            purrr_1.2.2            numDeriv_2016.8-1.1    jomo_2.7-6             e1071_1.7-17          
[21] mice_3.19.0            ggplot2_4.0.3          scales_1.4.0           iterators_1.0.14       MASS_7.3-65           
[26] cli_3.6.6              metafor_5.0-1          ragg_1.5.2             reformulas_0.4.4       generics_0.1.4        
[31] otel_0.2.0             rstudioapi_0.18.0      minqa_1.2.8            DBI_1.3.0              proxy_0.4-29          
[36] rlemon_0.2.1           stringr_1.6.0          operator.tools_1.6.3.1 splines_4.6.0          metadat_1.6-0         
[41] MetaUtility_2.1.2      mitools_2.4            vctrs_0.7.3            boot_1.3-32            glmnet_5.0            
[46] sandwich_3.1-1         hms_1.1.4              optmatch_0.10.8        mitml_0.4-5            systemfonts_1.3.2     
[51] foreach_1.5.2          tidyr_1.3.2            glue_1.8.1             nloptr_2.2.1           pan_2.0               
[56] chk_0.10.0             codetools_0.2-20       stringi_1.8.7          shape_1.4.6.1          gtable_0.3.6          
[61] lme4_2.0-1             tibble_3.3.1           pillar_1.11.1          R6_2.6.1               textshaping_1.0.5     
[66] Rdpack_2.6.6           arg_0.1.0              formula.tools_1.7.1    lattice_0.22-9         haven_2.5.5           
[71] rbibutils_2.4.1        backports_1.5.1        broom_1.0.13           class_7.3-23           Rcpp_1.1.1-1          
[76] nlme_3.1-169           mgcv_1.9-4             zoo_1.8-15             forcats_1.0.1          pkgconfig_2.0.3       
