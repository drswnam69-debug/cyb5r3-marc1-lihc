# =====================================================================
#  13_pancancer.R
#  목적 (심사 대응):
#   · R1-4 / R3 "Apply multiple-testing correction across the 27-cohort
#     pan-cancer scan, or state plainly that it is uncorrected."
#   · R3 "include the full 27-cohort … tables rather than deferring them to an
#     undeposited output directory"
#   · R3 "The pan-cancer LIHC estimate should not be called a 'reproduction' of
#     the raw-data Cox model, since the patients are the same."
#       → 출력 표에 overlap_with_primary 열을 두어 동일 환자임을 명시한다.
#
#  산출: results/13_pancancer_full.csv  (전 코호트, HR·CI·p·BH q, T-vs-N log2FC·p·q)
#        figures/13_pancancer_forest.png
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
hdr("13 · 범암종 주사 (다중검정 보정 포함)")

MIN_N <- 20; MIN_EVENTS <- 10; MIN_NORMAL <- 10

## 1. 발현 (PanCanAtlas EB++ 배치보정)
m <- fetch_gene(XH$pan, XDS$pan_expr, "CYB5R3", probemap = FALSE)
if (is.null(m)) stop("PanCanAtlas CYB5R3 조회 실패")
v <- tovec(m); names(v) <- nb(names(v))
EXP <- data.frame(sample = names(v), CYB5R3 = as.numeric(v))
EXP$patient <- pat(EXP$sample); EXP$stype <- typ(EXP$sample)

## 2. 임상/생존 (TCGA Clinical Data Resource) — 암종 약어 포함
##    ID 가 하드코딩이면 Xena 개정 시 깨지므로 XenaData 에서 찾아 쓴다.
sv <- xena_table("Survival_SupplementalTable", exact = XDS$pan_surv)
if (is.null(sv)) stop("TCGA-CDR 다운로드 실패 — 네트워크 또는 데이터셋 ID 확인")
message("  [diag] CDR 열이름: ", paste(utils::head(names(sv), 20), collapse = " | "))
gg <- function(p) { h <- grep(p, names(sv), ignore.case = TRUE, value = TRUE); if (length(h)) h[1] else NA }
col_pt <- gg("^sample$|bcr_patient_barcode|^_PATIENT$|^patient$")
col_ct <- gg("cancer.?type|^_primary_disease$|^type$|^disease$")
col_os <- gg("^OS$"); col_ot <- gg("^OS\\.?time$")
message(sprintf("  [diag] 매핑: patient='%s' cohort='%s' OS='%s' OS.time='%s'",
                col_pt, col_ct, col_os, col_ot))
if (any(is.na(c(col_pt, col_ct, col_os, col_ot))))
  stop("CDR 열 매핑 실패 — 위 [diag] 열이름 목록을 보내 주십시오.")
SV <- data.frame(patient = pat(sv[[col_pt]]),
                 cohort  = as.character(sv[[col_ct]]),
                 OS      = suppressWarnings(as.integer(sv[[col_os]])),
                 OS.time = suppressWarnings(as.numeric(sv[[col_ot]])))
SV <- SV[!is.na(SV$cohort) & nzchar(SV$cohort), ]
SV <- SV[!is.na(SV$OS) & !is.na(SV$OS.time) & SV$OS.time > 0, ]
SV <- SV %>% dplyr::group_by(patient) %>% dplyr::slice(1) %>% dplyr::ungroup()
message("  CDR patients: ", nrow(SV), " ; cohorts: ", dplyr::n_distinct(SV$cohort))

TUM <- EXP[EXP$stype %in% TUMOR_CODES, ] %>%
  dplyr::group_by(patient) %>% dplyr::slice(1) %>% dplyr::ungroup()
D <- dplyr::inner_join(TUM, SV, by = "patient")

## 3. 코호트별 단변량 Cox (per SD)
rows <- list()
for (co in sort(unique(D$cohort))) {
  d <- D[D$cohort == co & !is.na(D$CYB5R3), ]
  if (nrow(d) < MIN_N || sum(d$OS == 1) < MIN_EVENTS) next
  d$z <- as.numeric(scale(d$CYB5R3))
  s <- summary(coxph(Surv(OS.time, OS) ~ z, data = d))
  rows[[co]] <- data.frame(cohort = co, n = nrow(d), events = sum(d$OS == 1),
                           HR_perSD = s$conf.int[1,"exp(coef)"],
                           CI_low = s$conf.int[1,"lower .95"],
                           CI_high = s$conf.int[1,"upper .95"],
                           p_cox = s$coefficients[1,"Pr(>|z|)"])
}
COX <- dplyr::bind_rows(rows)
COX$q_cox <- p.adjust(COX$p_cox, method = "BH")     # ← 심사 요구: 코호트 전체 BH 보정

## 4. 코호트별 종양 vs 인접정상
NORM <- EXP[EXP$stype %in% NORMAL_CODES, ]
NORM <- dplyr::inner_join(NORM, SV[, c("patient","cohort")], by = "patient")
tn <- list()
for (co in sort(unique(D$cohort))) {
  tt <- D$CYB5R3[D$cohort == co]; nn <- NORM$CYB5R3[NORM$cohort == co]
  tt <- tt[!is.na(tt)]; nn <- nn[!is.na(nn)]
  if (length(nn) < MIN_NORMAL) next
  w <- wilcox.test(tt, nn)
  tn[[co]] <- data.frame(cohort = co, n_tumor = length(tt), n_normal = length(nn),
                         log2FC = median(tt) - median(nn), p_TvsN = w$p.value)
}
TN <- dplyr::bind_rows(tn)
if (nrow(TN)) TN$q_TvsN <- p.adjust(TN$p_TvsN, method = "BH")

OUT <- dplyr::left_join(COX, TN, by = "cohort")
OUT$direction_survival <- ifelse(OUT$HR_perSD > 1, "adverse", "protective")
OUT$significant_after_FDR <- OUT$q_cox < 0.05
OUT$overlap_with_primary <- ifelse(OUT$cohort == "LIHC",
  "SAME PATIENTS as the primary TCGA-LIHC analysis - not an independent replication",
  "")
OUT <- OUT[order(OUT$q_cox), ]
print(OUT, digits = 3); w_res(OUT, "13_pancancer_full.csv")

message(sprintf("\n  코호트 %d개 검정, BH 보정 후 유의 %d개.",
                nrow(OUT), sum(OUT$significant_after_FDR, na.rm = TRUE)))
if ("LIHC" %in% OUT$cohort) {
  l <- OUT[OUT$cohort == "LIHC", ]
  message(sprintf("  LIHC: HR %.2f (%.2f-%.2f), p=%.3g, BH q=%.3g  ← 주분석과 동일 환자",
                  l$HR_perSD, l$CI_low, l$CI_high, l$p_cox, l$q_cox))
}

try({
  pl <- OUT[!is.na(OUT$HR_perSD), ]
  pl$cohort <- factor(pl$cohort, levels = pl$cohort[order(pl$HR_perSD)])
  g <- ggplot2::ggplot(pl, ggplot2::aes(HR_perSD, cohort)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = CI_low, xmax = CI_high), height = .2) +
    ggplot2::geom_point(ggplot2::aes(colour = significant_after_FDR), size = 2) +
    ggplot2::scale_colour_manual(values = c("TRUE" = "#C0392B", "FALSE" = "grey55"),
                                 name = "BH q < 0.05") +
    ggplot2::scale_x_log10() + ggplot2::theme_minimal() +
    ggplot2::labs(x = "Overall-survival HR per SD of CYB5R3 (univariable)", y = NULL,
                  title = "CYB5R3 prognostic direction across TCGA cohorts (BH-corrected)")
  ggplot2::ggsave(file.path(FIG, "13_pancancer_forest.png"), g, width = 7, height = 8, dpi = 300)
}, silent = TRUE)

save_session("13_pancancer")
message("\n[13] 완료.")
