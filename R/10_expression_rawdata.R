# =====================================================================
#  10_expression_rawdata.R
#  목적 (심사 대응):
#   · R1 minor 2, R3 "Rebuild the primary analyses from raw data"
#     → 원 투고 §2.1–2.5 / Figure 2 / Table 1 의 '그림에서 눈으로 읽은' 발현값을
#       전부 원자료에서 재생성한다. 정확한 n, 중앙값, IQR, log2FC, Wilcoxon p,
#       Cliff's delta 를 산출한다.
#   · R3 "Reconcile the figures that contradict one another"
#     → 동일 질문(LIHC 종양 vs 정상)에 대해 서로 다른 세 추정치(UALCAN /
#       Xena TOIL / PanCanAtlas EB++)를 '대조군·전처리'와 함께 한 표에 나란히
#       제시한다. 어느 것이 primary 인지 표에 명시된다.
#
#  산출: results/10_expression_estimates.csv
#        results/10_expression_reconciliation.csv
#        figures/10_tumor_vs_normal_rawdata.png
# =====================================================================
.here <- function() {
  fr <- sys.frames()
  if (length(fr)) for (i in rev(seq_along(fr))) {
    o <- fr[[i]]$ofile
    if (!is.null(o) && nzchar(o)) return(dirname(normalizePath(o)))
  }
  m <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(m)) return(dirname(normalizePath(sub("^--file=", "", m[1]))))
  getwd()
}
if (!exists("RES")) source(file.path(.here(), "01_common.R"))
hdr("10 · 발현 원자료 재생성")

GENES <- AXIS
out   <- list()

## =====================================================================
## (A) TCGA HiSeqV2 (log2(norm_count+1)) — 종양 vs 인접정상, LIHC & LUAD
##     주의: MTARC1 은 legacy probeMap 에 MARC1/MOSC1 로 색인될 수 있음
## =====================================================================
hdr("A. TCGA HiSeqV2 — 종양 vs 인접정상")
for (co in c("LIHC","LUAD")) {
  ds <- XDS[[paste0(tolower(co), "_hiseq")]]
  for (g in GENES) {
    m <- fetch_gene(XH$tcga, ds, g)
    if (is.null(m)) { message("  [skip] ", co, "/", g); next }
    v <- tovec(m); names(v) <- nb(names(v))
    d <- data.frame(sample = names(v), val = as.numeric(v))
    d$grp <- ifelse(typ(d$sample) %in% TUMOR_CODES,  "Tumor",
             ifelse(typ(d$sample) %in% NORMAL_CODES, "Normal", NA))
    d <- d[!is.na(d$grp) & !is.na(d$val), ]
    tt <- d$val[d$grp == "Tumor"]; nn <- d$val[d$grp == "Normal"]
    if (length(nn) < 5) { message("  [skip] ", co, "/", g, ": 정상 표본 부족"); next }
    w <- wilcox.test(tt, nn)
    out[[length(out)+1]] <- data.frame(
      cohort = co, gene = g,
      dataset = "TCGA HiSeqV2 (log2 norm_count+1)",
      normal_reference = "TCGA adjacent normal (sample type 11)",
      alias_used = attr(m, "alias_used"),
      n_tumor = length(tt), n_normal = length(nn),
      med_tumor = median(tt), q1_tumor = quantile(tt,.25), q3_tumor = quantile(tt,.75),
      med_normal = median(nn), q1_normal = quantile(nn,.25), q3_normal = quantile(nn,.75),
      log2FC_median = median(tt) - median(nn),
      cliffs_delta = cliffs_delta(tt, nn),
      wilcox_p = w$p.value)
    message(sprintf("  %s/%s  T n=%d med=%.3f | N n=%d med=%.3f | log2FC=%+.3f p=%.3g",
                    co, g, length(tt), median(tt), length(nn), median(nn),
                    median(tt)-median(nn), w$p.value))
  }
}

## =====================================================================
## (B) Xena TOIL TcgaTargetGtex (log2(TPM+0.001)) — 종양 vs TCGA인접정상 vs GTEx정상
##     이 데이터셋은 TCGA 와 GTEx 를 동일 파이프라인으로 재처리한 것이므로
##     '조직 간 / 코호트 간' 비교에 적합하다. 선형 TPM 중앙값도 함께 보고한다.
## =====================================================================
hdr("B. Xena TOIL — 종양 vs 인접정상 vs GTEx 정상")
ph_file <- file.path(CACHE, "TcgaTargetGTEX_phenotype.txt.gz")
if (!file.exists(ph_file))
  utils::download.file(sprintf("%s/download/%s.gz", XH$toil, XDS$toil_pheno),
                       ph_file, mode = "wb", quiet = TRUE)
ph <- data.table::fread(ph_file)
data.table::setnames(ph, 1, "sample")
gcol <- function(p) grep(p, names(ph), ignore.case = TRUE, value = TRUE)[1]
ph2 <- data.frame(
  sample = nb(ph$sample),
  site   = as.character(ph[[gcol("_primary_site|primary.site")]]),
  stype  = as.character(ph[[gcol("_sample_type|sample.type")]]),
  study  = as.character(ph[[gcol("_study|^study$")]]),
  detail = if (!is.na(gcol("detailed|primary.disease|disease.or.tissue")))
             as.character(ph[[gcol("detailed|primary.disease|disease.or.tissue")]]) else NA)

toil_group <- function(d, tissue, disease_pat) {
  dplyr::case_when(
    d$study == "GTEX" & grepl(tissue, d$site, ignore.case = TRUE) ~ "GTEx normal",
    d$study == "TCGA" & grepl(disease_pat, d$detail, ignore.case = TRUE) &
      grepl("Primary Tumor", d$stype, ignore.case = TRUE)         ~ "TCGA tumor",
    d$study == "TCGA" & grepl(disease_pat, d$detail, ignore.case = TRUE) &
      grepl("Normal Tissue|Solid Tissue Normal", d$stype, ignore.case = TRUE) ~ "TCGA adjacent normal",
    TRUE ~ NA_character_)
}
TOIL_DEF <- list(
  LIHC = list(tissue = "Liver", disease = "Hepatocellular"),
  LUAD = list(tissue = "Lung",  disease = "Lung Adenocarcinoma")
)
for (co in names(TOIL_DEF)) {
  for (g in GENES) {
    m <- fetch_gene(XH$toil, XDS$toil_tpm, g)
    if (is.null(m)) { message("  [skip] TOIL/", g); next }
    v <- tovec(m); names(v) <- nb(names(v))
    D <- merge(data.frame(sample = names(v), val = as.numeric(v)), ph2, by = "sample")
    D$grp <- toil_group(D, TOIL_DEF[[co]]$tissue, TOIL_DEF[[co]]$disease)
    D <- D[!is.na(D$grp) & !is.na(D$val), ]
    tt <- D$val[D$grp == "TCGA tumor"]
    for (ref in c("TCGA adjacent normal","GTEx normal")) {
      nn <- D$val[D$grp == ref]
      if (length(nn) < 5 || length(tt) < 5) next
      w <- wilcox.test(tt, nn)
      out[[length(out)+1]] <- data.frame(
        cohort = co, gene = g,
        dataset = "Xena TOIL TcgaTargetGtex RSEM (log2 TPM+0.001)",
        normal_reference = ref,
        alias_used = attr(m, "alias_used"),
        n_tumor = length(tt), n_normal = length(nn),
        med_tumor = median(tt), q1_tumor = quantile(tt,.25), q3_tumor = quantile(tt,.75),
        med_normal = median(nn), q1_normal = quantile(nn,.25), q3_normal = quantile(nn,.75),
        log2FC_median = median(tt) - median(nn),
        cliffs_delta = cliffs_delta(tt, nn),
        wilcox_p = w$p.value)
      message(sprintf("  TOIL %s/%s vs %-22s T n=%d med=%.3f | N n=%d med=%.3f | log2FC=%+.3f p=%.3g",
                      co, g, ref, length(tt), median(tt), length(nn), median(nn),
                      median(tt)-median(nn), w$p.value))
    }
    ## 원고 본문용 선형 TPM 중앙값 (그림에서 읽은 '~120→170 TPM' 을 대체)
    lin <- function(x) 2^x - 0.001
    message(sprintf("     [linear TPM] tumor median=%.1f | adjacent normal median=%.1f | GTEx median=%.1f",
                    median(lin(tt)),
                    suppressWarnings(median(lin(D$val[D$grp=="TCGA adjacent normal"]))),
                    suppressWarnings(median(lin(D$val[D$grp=="GTEx normal"])))))
  }
}

## =====================================================================
## (C) PanCanAtlas EB++ batch-corrected — 보충 Table S3 이 쓴 자료
##     (원 투고에서 LIHC log2FC −0.07, p=0.088 이 나온 출처)
## =====================================================================
hdr("C. PanCanAtlas EB++ (Supplementary Table S3 출처)")
try({
  sv_pheno <- NULL
  for (co in c("LIHC","LUAD")) for (g in GENES) {
    m <- fetch_gene(XH$pan, XDS$pan_expr, g, probemap = FALSE)
    if (is.null(m)) { message("  [skip] PanCan/", g); next }
    v <- tovec(m); names(v) <- nb(names(v))
    ## PanCanAtlas 발현은 전체 암종을 포함 → 코호트 배정을 위해 TOIL 표현형 사용
    D <- merge(data.frame(sample = bc15(names(v)), val = as.numeric(v)),
               data.frame(sample = bc15(ph2$sample), detail = ph2$detail,
                          stype = ph2$stype, study = ph2$study),
               by = "sample")
    D <- D[D$study == "TCGA", ]
    pat_d <- if (co == "LIHC") "Hepatocellular" else "Lung Adenocarcinoma"
    D <- D[grepl(pat_d, D$detail, ignore.case = TRUE), ]
    D$grp <- ifelse(grepl("Primary Tumor", D$stype, ignore.case = TRUE), "Tumor",
             ifelse(grepl("Normal", D$stype, ignore.case = TRUE), "Normal", NA))
    D <- D[!is.na(D$grp) & !is.na(D$val), ]
    tt <- D$val[D$grp == "Tumor"]; nn <- D$val[D$grp == "Normal"]
    if (length(nn) < 5) next
    w <- wilcox.test(tt, nn)
    out[[length(out)+1]] <- data.frame(
      cohort = co, gene = g,
      dataset = "Xena PanCanAtlas EB++ batch-corrected (log2)",
      normal_reference = "TCGA adjacent normal (same cohort)",
      alias_used = attr(m, "alias_used"),
      n_tumor = length(tt), n_normal = length(nn),
      med_tumor = median(tt), q1_tumor = quantile(tt,.25), q3_tumor = quantile(tt,.75),
      med_normal = median(nn), q1_normal = quantile(nn,.25), q3_normal = quantile(nn,.75),
      log2FC_median = median(tt) - median(nn),
      cliffs_delta = cliffs_delta(tt, nn),
      wilcox_p = w$p.value)
    message(sprintf("  PanCan %s/%s  T n=%d | N n=%d | log2FC=%+.3f p=%.3g",
                    co, g, length(tt), length(nn), median(tt)-median(nn), w$p.value))
  }
}, silent = FALSE)

## =====================================================================
## 저장 + 화해(reconciliation) 표
## =====================================================================
res <- dplyr::bind_rows(out)
if (nrow(res)) {
  res <- res %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, 6)))
  w_res(res, "10_expression_estimates.csv")

  rec <- res %>%
    dplyr::filter(cohort == "LIHC", gene == "CYB5R3") %>%
    dplyr::transmute(
      estimate = paste(dataset, "vs", normal_reference),
      n = sprintf("%d tumor / %d normal", n_tumor, n_normal),
      log2FC = sprintf("%+.3f", log2FC_median),
      p = signif(wilcox_p, 3),
      cliffs_delta = round(cliffs_delta, 3))
  rec$note <- "UALCAN (TCGA level-3 TPM, unpaired t-test) is a portal-rendered value and is entered as a separate row in Table 2; it is not computed by this script."
  w_res(rec, "10_expression_reconciliation.csv")
  print(rec)
} else message("[10] 결과 없음 — 위 진단 메시지를 확인하십시오.")

save_session("10_expression")
message("\n[10] 완료.")
