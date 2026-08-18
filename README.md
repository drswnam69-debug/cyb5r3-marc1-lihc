# CYB5R3–mARC1 redox axis in MASLD and hepatocellular carcinoma — analysis code

Reproducible analysis pipeline accompanying the manuscript

> Nam SW. *The CYB5R3–mARC1 redox axis across the human MASLD-to-hepatocellular-carcinoma
> continuum.* Submitted to *International Journal of Molecular Sciences* (ijms-4481300),
> revision 1.

All analyses use publicly available, de-identified human data. No new human or animal
data were generated.

---

## 1. What this repository does

| Script | Question it answers | Key outputs |
|---|---|---|
| `R/01_common.R` | Paths, packages, statistical helpers (shared; not run directly) | — |
| `R/10_expression_rawdata.R` | Tumour-versus-normal expression regenerated **from raw data**, with three independent estimates placed side by side | `10_expression_estimates.csv`, `10_expression_reconciliation.csv` |
| `R/11_survival_primary.R` | **Primary prognostic analysis**: continuous per-SD multivariable Cox, patient flow, events, follow-up, proportional-hazards diagnostics, non-linearity, C-index and ΔC-index, joint/interaction models, extended covariates, multiple imputation, covariate-adjusted external replication | `11_flow.csv`, `11_cox_primary.csv`, `11_cox_ph.csv`, `11_cindex.csv`, `11_cox_joint.csv`, `11_gse14520_adjusted.csv` … |
| `R/12_coexpression_tissue.R` | Full 16-gene panel with 95% CIs and q-values; **formal Fisher r-to-z test** of the liver-versus-lung difference; partial correlations adjusted for tumour cellular composition; permutation null for the number of sign reversals | `12_coexpression_full.csv`, `12_tissue_difference.csv`, `12_purity_adjusted.csv`, `12_permutation_reversals.csv` |
| `R/13_pancancer.R` | Pan-cancer scan across all TCGA cohorts **with Benjamini–Hochberg correction**; flags that the LIHC estimate uses the same patients as the primary analysis | `13_pancancer_full.csv` |
| `R/14_masld_cohorts.R` | The axis across **three human MASLD liver-biopsy cohorts** with pathologist-assigned histology (GSE130970, GSE135251, GSE162694), FDR-corrected, with random-effects meta-analysis | `14_masld_trend_all.csv`, `14_masld_meta.csv`, `14_masld_partial.csv` |
| `R/15_depmap.R` | CRISPR dependency of CYB5R3 in liver versus lung cell lines (a **negative** result, reported as such) | `15_depmap_dependency.csv` |
| `R/16_lihc_stratified.R` | Stage, histologic grade and viral-serology stratification within TCGA-LIHC | `16_lihc_stage_grade.csv`, `16_lihc_etiology.csv`, `16_lihc_serologynegative_cox.csv` |
| `R/17_manuscript_figures.R` | Publication figures at 300 dpi (PNG + LZW TIFF), each with its numerical source data | `figures/Figure*.png/.tiff`, `figures/Figure*_source.csv` |
| `R/run_all.R` | Runs everything in order and writes a timestamped log | `run_all_status.csv` |

For convenience, three one-click files sit at the top level: **`RUN_ALL.R`** (whole pipeline),
**`RUN_FIGURES.R`** (figures only) and **`RUN_16.R`** (one step). Open one in RStudio and press Source.
**`RESULTS_SUMMARY.md`** lists every headline number with the file it came from.

Every script writes `sessionInfo_<tag>.txt` and `RUNLOG_<tag>.txt` to `results/`.
The random seed is fixed at `20260817` in `R/01_common.R`.

---

## 2. Data sources

| Resource | Dataset identifier | Access |
|---|---|---|
| UCSC Xena — TCGA hub | `TCGA.LIHC.sampleMap/HiSeqV2`, `TCGA.LUAD.sampleMap/HiSeqV2`, `TCGA.LIHC.sampleMap/LIHC_clinicalMatrix` | `UCSCXenaTools` |
| UCSC Xena — TOIL hub | `TcgaTargetGtex_rsem_gene_tpm`, `TcgaTargetGTEX_phenotype.txt` | `UCSCXenaTools` |
| UCSC Xena — PanCanAtlas hub | `EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena`, `Survival_SupplementalTable_S1_20171025_xena_sp`, `TCGA_mastercalls.abs_tables_JSedit.fixed.txt` (ABSOLUTE tumour purity) | `UCSCXenaTools` |
| GEO | GSE14520 (HCC, external survival replication); GSE130970, GSE135251, GSE162694 (human MASLD liver biopsies) | `GEOquery` + NCBI FTP |
| DepMap | CRISPR (Chronos) gene effect | Bioconductor `depmap` / ExperimentHub, or portal CSV |

Primary publications for the human cohorts:
Roessler S, et al. *Cancer Res* 2010;70:10202–12 (GSE14520);
Hoang SA, et al. *Sci Rep* 2019;9:12541 (GSE130970);
Govaere O, et al. *Sci Transl Med* 2020;12:eaba4448 (GSE135251);
Pantano L, et al. *Sci Rep* 2021;11:18045 (GSE162694).

---

## 3. 실행 방법 (교수님용)

### 3.1 준비

```r
# R 4.3 이상 권장
install.packages(c("UCSCXenaTools","data.table","dplyr","tidyr","stringr","readr",
                   "tibble","purrr","ggplot2","ggpubr","survival","survminer",
                   "boot","R.utils","rms","mice","metafor"),
                 repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery","Biobase","depmap","ExperimentHub"), ask = FALSE, update = FALSE)
```

`rms`, `mice`, `metafor` 는 **선택**입니다. 없으면 해당 분석만 건너뛰고 나머지는 정상 실행됩니다.
다만 `mice`(다중대치)와 `metafor`(메타분석)는 심사 응답에 쓰이므로 설치를 권합니다.

### 3.2 실행

RStudio에서 `R/run_all.R` 을 열고 **Source** 하시거나, 터미널에서:

```bash
cd <이 폴더>
Rscript R/run_all.R
```

한 단계만 다시 돌리려면:

```r
RUN_ONLY <- c("14_masld_cohorts.R")
source("R/run_all.R")
```

### 3.3 이미 받아 둔 GEO 자료 재사용

`Revision_2026_human_MASLD/analysis/cache/geo/` 에 GSE130970·GSE135251·GSE162694 를
이미 내려받아 두셨습니다. 스크립트가 자동으로 그 경로를 찾습니다. 경로가 다르면:

```r
Sys.setenv(CYB5R3_GEO_CACHE = "/전체/경로/cache/geo")
source("R/run_all.R")
```

### 3.4 예상 소요시간

| 단계 | 대략 |
|---|---|
| 10 · 발현 재생성 | 5–15분 (TOIL 대용량 슬라이스) |
| 11 · 생존분석 | 5–10분 (부트스트랩 1,000회 × 3) |
| 12 · 공발현 | 10–20분 (순열 10,000회) |
| 13 · 범암종 | 5–10분 |
| 14 · MASLD | 캐시 있으면 3–5분, 새로 받으면 20분+ |
| 15 · DepMap | 2–5분 (ExperimentHub 최초 1회 다운로드 시 더 걸림) |

### 3.5 실행이 끝나면

`results/` 폴더 전체와 `figures/` 폴더를 저에게 보내 주시면, 그 수치로 본문·표·응답서를
확정하겠습니다. 특히 아래 파일이 원고의 핵심 수치를 담습니다.

- `11_flow.csv` — 환자 흐름·사망수·추적기간 (심사 R3 필수 요구)
- `11_cox_primary.csv`, `11_cindex.csv` — 주분석
- `12_tissue_difference.csv`, `12_permutation_reversals.csv` — 조직 역전의 형식적 검정
- `14_masld_meta.csv` — MASLD 통합 결과
- `10_expression_reconciliation.csv` — 세 추정치 화해표

---

## 4. Reproducibility notes

* A single random seed (`20260817`) is set in `R/01_common.R` and reported in every run log.
* Every result file is written by exactly one script; no numbers in the manuscript are
  transcribed from rendered portal plots.
* Where a portal value is quoted (UALCAN, KM Plotter, cBioPortal, CPTAC, Human Protein
  Atlas), it is labelled as a portal-derived observation in the manuscript and is **not**
  presented as an independent analysis. UALCAN, cBioPortal, KM Plotter's liver cohort and
  the pan-cancer scan all draw on TCGA; GSE14520 is the only independent survival cohort.
* Session information and package versions for each run are written to `results/`.

## 5. Licence

Code: MIT (see `LICENSE`). The underlying data remain subject to the terms of their
respective repositories (TCGA/GDC, GTEx, CPTAC, GEO, DepMap).

## 6. Citation

If you use this code, please cite the manuscript and this archive (Zenodo DOI printed on
release).
