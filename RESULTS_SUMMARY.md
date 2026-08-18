# Results summary

Headline numbers produced by this pipeline, with the file each comes from. Every value
in the manuscript is taken from these files; nothing is transcribed from a rendered plot.

Run: 2026-08-18. R 4.6.1 (aarch64-apple-darwin23). Seed 20260817. See `results/sessionInfo_*.txt`.

## Primary prognostic analysis — TCGA-LIHC (`results/11_*.csv`)

| Quantity | Value |
|---|---|
| Patient flow | 423 samples → 371 primary tumours → 369 matched to *MTARC1* → 363 with survival → **339 analysed** |
| Deaths | **114** |
| Median follow-up (reverse Kaplan–Meier) | **816 days (26.8 months)** |
| *CYB5R3*, per SD, adjusted for age/sex/stage | **HR 1.226 (1.015–1.480), p = 0.035** |
| *MTARC1*, per SD, adjusted | HR 0.936 (0.798–1.098), p = 0.420 |
| Progression-free interval (162 events) | *CYB5R3* 1.007, p = 0.93; *MTARC1* 1.017, p = 0.82 |
| Proportional hazards | *CYB5R3* term p = 0.91; global p = 0.34 |
| Non-linearity (spline vs linear) | LR p = 0.51; AIC 1134.8 vs 1137.5 |
| **C-index, clinical base model** | **0.629 (0.580–0.688)** |
| **ΔC-index on adding *CYB5R3*** | **0.005 (−0.015 to +0.037)** |
| Multiple imputation (m = 20, n = 363, 128 deaths) | HR 1.232 (1.026–1.479), p = 0.026 |
| Joint model, *MTARC1* added to *CYB5R3* | LR χ² = 0.75, 1 df, **p = 0.39** |
| Interaction *CYB5R3* × *MTARC1* | HR 1.039, p = 0.57 |
| **Sample-level correlation of the two transcripts** | **rₛ = −0.027 (−0.133 to +0.080), p = 0.62** |
| + log10(AFP) (n = 261, 74 deaths) | HR 1.165 (0.923–1.471), p = 0.199 |
| + vascular invasion (n = 291, 85 deaths) | HR 1.156 (0.926–1.442), p = 0.201 |
| + cirrhosis, Ishak 5–6 (n = 196, 58 deaths) | HR 1.171 (0.887–1.545), p = 0.266 |

## External cohort — GSE14520 (`results/11_gse14520_adjusted.csv`)

n = 242, 96 deaths (unadjusted); n = 225, 86 deaths (adjusted).

| Model | *CYB5R3* | *MTARC1* |
|---|---|---|
| Unadjusted | 1.420 (1.144–1.764), p = 0.0015 | 0.728 (0.597–0.888), p = 0.0017 |
| + age, sex, TNM | 1.244 (0.992–1.559), **p = 0.058** | 0.894 (0.725–1.102), **p = 0.294** |
| + cirrhosis, multinodularity | 1.218 (0.970–1.529), p = 0.089 | 0.881 (0.713–1.087), p = 0.237 |
| Joint (adjusted) | **1.267 (1.012–1.586), p = 0.039** | 0.865 (0.699–1.071), p = 0.184 |

Adjustment removes conventional significance from both single-gene estimates. This cohort
replicates the **direction** of effect, not independent prognostic value.

## Tumour-versus-normal expression (`results/10_expression_estimates.csv`)

The direction in liver depends on the processing pipeline, so no claim of over-expression is made.

| Cohort | Matrix | Normal reference | log2FC | p |
|---|---|---|---|---|
| LIHC | Xena TOIL | GTEx normal liver | +0.213 | 0.0016 |
| LIHC | Xena TOIL | TCGA adjacent normal | +0.480 | 7.1 × 10⁻⁷ |
| LIHC | TCGA HiSeqV2 | TCGA adjacent normal | −0.070 | 0.088 |
| LIHC | PanCanAtlas EB++ | TCGA adjacent normal | −0.060 | 0.102 |
| LUAD | all four | all references | −0.73 to −0.93 | 2.6 × 10⁻²⁴ to 1.2 × 10⁻⁶³ |

## Co-expression, liver versus lung (`results/12_*.csv`)

* Within LIHC, only *SCD* survives FDR correction (rₛ = +0.164, q = 0.024).
* Fisher r-to-z: six genes differ significantly between tissues after correction —
  *CAT* (q = 1.7 × 10⁻⁵), *PARP16* (0.0024), *TXNRD1* (0.0069), *SOD2* (0.0070),
  *SCD* (0.0081), *KEAP1* (0.021). **All six survive adjustment for cellular composition.**
* Permutation null: 8 sign reversals observed, 8.0 expected, **p = 0.58**. The *count* of
  reversals carries no information; only the specific coefficient differences do.
* ABSOLUTE consensus purity could not be retrieved from the Xena hubs; a stromal and immune
  marker principal component was used as a composition proxy, and the manuscript says so.

## Pan-cancer scan (`results/13_pancancer_full.csv`)

28 cohorts. **After Benjamini–Hochberg correction across cohorts, no cohort is significant**,
including LIHC (nominal HR 1.204, p = 0.044, **q = 0.194**). Smallest corrected value q = 0.146.

## Human MASLD liver biopsies (`results/14_masld_*.csv`)

GSE130970 n = 78, GSE135251 n = 216, GSE162694 n = 99–102 matched. Random-effects pooling,
Benjamini–Hochberg correction across all 72 gene × histology × cohort tests.

| Gene | Feature | GSE130970 | GSE135251 | GSE162694 | Pooled (95% CI) | q | I² |
|---|---|---|---|---|---|---|---|
| *CYB5R3* | fibrosis | +0.246 | +0.010 | **−0.221** | +0.008 (−0.243 to +0.258) | 0.95 | **83%** |
| *CYB5R3* | NAS | +0.333 | +0.127 | −0.131 | +0.110 (−0.147 to +0.353) | 0.50 | **83%** |
| ***MTARC1*** | **fibrosis** | −0.323 | −0.180 | −0.243 | **−0.224 (−0.317 to −0.128)** | **3 × 10⁻⁵** | **0%** |
| *PARP16* | fibrosis | −0.012 | −0.187 | −0.246 | −0.169 (−0.264 to −0.070) | 0.0023 | 0% |
| *SCD* | fibrosis | +0.016 | +0.138 | +0.202 | +0.131 (+0.032 to +0.228) | 0.023 | 0% |
| *NQO1* | fibrosis | +0.120 | +0.446 | +0.013 | +0.211 (−0.073 to +0.463) | 0.21 | 86% |
| *COL1A1* | fibrosis | +0.425 | +0.460 | +0.212 | +0.377 (+0.220 to +0.515) | 3 × 10⁻⁵ | 62% |

**The central observation.** *CYB5R3* estimates differ in sign between cohorts (I² = 83%),
while *MTARC1*, *PARP16* and *SCD*, measured in the same biopsies against the same
pathologist-assigned scores, each give a consistent pooled estimate with I² = 0%. The
inconsistency is therefore specific to *CYB5R3* rather than a property of the cohorts.

## Stage, grade and aetiology (`results/16_*.csv`)

* AJCC stage trend: Jonckheere–Terpstra **p = 0.164** (I 171, II 86, III 85, IV 5).
* Histologic grade: **p = 0.0067**; expression *declines* with worsening grade
  (G1 12.84, G2 12.85, G3 12.71, G4 12.37).
* **Aetiology.** The TCGA-LIHC `viral_hepatitis_serology` field lists the serologies that
  returned a **positive** result. It is populated for 164 of 371 tumours; an empty field
  means no positive serology was recorded and cannot be distinguished from "not tested".
  **No metabolic subset can be defined from these data.** Expression does not differ between
  the groups (p = 0.209). Exploratory survival: no positive serology recorded
  HR 1.297 (0.985–1.709), p = 0.064, n = 203/59 deaths; viral serology positive
  HR 1.130 (0.880–1.452), p = 0.337, n = 162/71 deaths.

## DepMap (`results/15_depmap_dependency.csv`)

Liver 22 lines, median Chronos −0.046, 0% below −0.5. Lung 133 lines, median −0.018,
0.75% below −0.5. **Wilcoxon p = 0.165.** *CYB5R3* is not a genetic dependency in either lineage.
