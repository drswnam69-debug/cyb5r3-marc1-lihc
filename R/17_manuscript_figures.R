# =====================================================================
#  17_manuscript_figures.R
#  원고용 그림 생성 (MDPI 규격: 폭 180 mm, 300 dpi, PNG + TIFF/LZW)
#
#  심사 대응:
#   · R1 minor 2 "Figure 2 should use raw-data distributions and exact
#     summary statistics instead of approximate median values read from
#     rendered plots"  → Figure 2 를 원자료 분포로 다시 그린다.
#   · R3 "numerical source data for every figure"  → 모든 그림은
#     results/*.csv 또는 원자료에서 직접 생성되며, 각 그림의 수치는
#     figures/<그림>_source.csv 로 함께 저장된다.
#   · R3 "Remove 'rises along the continuum'"  → Figure 5 는 코호트별
#     추정치와 I² 를 함께 보여 이질성을 감추지 않는다.
#
#  산출: figures/Figure2_expression.{png,tiff} + _source.csv
#        figures/Figure3B_coexpression.{png,tiff} + _source.csv
#        figures/Figure5_masld.{png,tiff} + _source.csv
#        figures/Figure6_prognostic.{png,tiff} + _source.csv
#        figures/FigureS6_patient_flow.{png,tiff}
#        figures/FigureS7_schoenfeld.{png,tiff}
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
hdr("17 · 원고용 그림 생성")

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

## ---- 검증된 색상 슬롯 (색각이상 안전; 별도 검증 통과) ------------------
PAL <- c(blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a",
         yellow = "#eda100", violet = "#4a3aa7", grey = "#8a8a85")
INK <- c(primary = "#0b0b0b", secondary = "#52514e", muted = "#8a8a85")

theme_ms <- function(base = 9) {
  theme_minimal(base_size = base, base_family = "") +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
      axis.text          = element_text(colour = INK["secondary"]),
      axis.title         = element_text(colour = INK["primary"]),
      strip.text         = element_text(colour = INK["primary"], face = "bold", size = base),
      plot.title         = element_text(colour = INK["primary"], face = "bold", size = base + 2),
      plot.subtitle      = element_text(colour = INK["secondary"], size = base - 0.5),
      plot.caption       = element_text(colour = INK["muted"], size = base - 1.5, hjust = 0),
      legend.position    = "top",
      legend.title       = element_blank(),
      legend.text        = element_text(colour = INK["secondary"]),
      plot.background    = element_rect(fill = "white", colour = NA),
      panel.background   = element_rect(fill = "white", colour = NA))
}

MM <- function(mm) mm / 25.4          # mm → inch
save_fig <- function(p, name, w_mm = 180, h_mm = 110) {
  png_p  <- file.path(FIG, paste0(name, ".png"))
  tiff_p <- file.path(FIG, paste0(name, ".tiff"))
  ggsave(png_p,  p, width = MM(w_mm), height = MM(h_mm), dpi = 300, units = "in")
  ggsave(tiff_p, p, width = MM(w_mm), height = MM(h_mm), dpi = 300, units = "in",
         compression = "lzw", device = "tiff")
  message("  → ", png_p)
  message("  → ", tiff_p)
}
save_src <- function(df, name) {
  p <- file.path(FIG, paste0(name, "_source.csv"))
  data.table::fwrite(df, p); message("  → ", p)
}
rd <- function(f) as.data.frame(data.table::fread(file.path(RES, f)))
fmt_p <- function(p) ifelse(is.na(p), "NA",
                     ifelse(p < 1e-4, sprintf("%.0e", p), sprintf("%.3f", p)))

## =====================================================================
## Figure 2 — 종양 vs 정상 발현 (원자료 분포)
## =====================================================================
hdr("Figure 2 · 발현 원자료 분포")
try({
  ph_file <- file.path(CACHE, "TcgaTargetGTEX_phenotype.txt.gz")
  if (!file.exists(ph_file))
    utils::download.file(sprintf("%s/download/%s.gz", XH$toil, XDS$toil_pheno),
                         ph_file, mode = "wb", quiet = TRUE)
  ph <- data.table::fread(ph_file); data.table::setnames(ph, 1, "sample")
  gcol <- function(p) grep(p, names(ph), ignore.case = TRUE, value = TRUE)[1]
  P <- data.frame(sample = nb(ph$sample),
                  site   = as.character(ph[[gcol("_primary_site|primary.site")]]),
                  stype  = as.character(ph[[gcol("_sample_type|sample.type")]]),
                  study  = as.character(ph[[gcol("_study|^study$")]]),
                  detail = as.character(ph[[gcol("detailed|primary.disease|disease.or.tissue")]]))

  DEF <- list(LIHC = list(tissue = "Liver", disease = "Hepatocellular"),
              LUAD = list(tissue = "Lung",  disease = "Lung Adenocarcinoma"))
  rows <- list()
  for (g in AXIS) {
    m <- fetch_gene(XH$toil, XDS$toil_tpm, g)
    if (is.null(m)) { message("  [skip] ", g); next }
    v <- tovec(m); names(v) <- nb(names(v))
    D0 <- merge(data.frame(sample = names(v), val = as.numeric(v)), P, by = "sample")
    for (co in names(DEF)) {
      D <- D0
      D$grp <- dplyr::case_when(
        D$study == "GTEX" & grepl(DEF[[co]]$tissue, D$site, ignore.case = TRUE) ~ "GTEx normal",
        D$study == "TCGA" & grepl(DEF[[co]]$disease, D$detail, ignore.case = TRUE) &
          grepl("Primary Tumor", D$stype, ignore.case = TRUE) ~ "Tumor",
        D$study == "TCGA" & grepl(DEF[[co]]$disease, D$detail, ignore.case = TRUE) &
          grepl("Normal", D$stype, ignore.case = TRUE) ~ "Adjacent normal",
        TRUE ~ NA_character_)
      D <- D[!is.na(D$grp) & !is.na(D$val), c("sample","val","grp")]
      if (nrow(D)) rows[[length(rows)+1]] <- cbind(cohort = co, gene = g, D)
    }
  }
  E <- dplyr::bind_rows(rows)
  E$grp <- factor(E$grp, levels = c("GTEx normal","Adjacent normal","Tumor"))
  E$panel <- factor(paste0(E$gene, "  —  ", E$cohort),
                    levels = c("CYB5R3  —  LIHC","CYB5R3  —  LUAD",
                               "MTARC1  —  LIHC","MTARC1  —  LUAD"))

  ## 각 패널의 검정 결과 (종양 vs 각 정상군)
  ann <- E %>% dplyr::group_by(panel, cohort, gene) %>% dplyr::group_modify(function(d, k) {
    tt <- d$val[d$grp == "Tumor"]; out <- list()
    for (ref in c("GTEx normal","Adjacent normal")) {
      nn <- d$val[d$grp == ref]
      if (length(nn) < 5 || length(tt) < 5) next
      w <- suppressWarnings(wilcox.test(tt, nn))
      out[[ref]] <- data.frame(reference = ref, n_tumor = length(tt), n_normal = length(nn),
                               log2FC = median(tt) - median(nn), p = w$p.value)
    }
    dplyr::bind_rows(out)
  }) %>% dplyr::ungroup()

  lab <- ann %>% dplyr::group_by(panel) %>% dplyr::summarise(
    txt = paste(sprintf("vs %s: %+.2f, p = %s", sub(" normal","", reference), log2FC, fmt_p(p)),
                collapse = "\n"), .groups = "drop")
  nlab <- E %>% dplyr::count(panel, grp)

  ## ---- 축 압축 문제 해결 -------------------------------------------
  ## LUAD 패널에는 TPM = 0 인 소수의 검체가 있어 log2(TPM+0.001) 이 약 -10 이
  ## 되고, 그 결과 상자그림이 패널 높이의 10 % 남짓으로 눌린다.
  ## 상자·수염 통계량은 '모든' 자료로 계산하되(따라서 요약값은 바뀌지 않는다),
  ## 표시 범위는 수염 바깥으로만 잡고 범위를 벗어난 개별 점은 그리지 않는다.
  ## 그리지 않은 점의 개수는 그림 각주와 _source.csv 에 함께 보고한다.
  bx <- E %>% dplyr::group_by(panel, grp) %>% dplyr::group_modify(function(d, k) {
    q <- stats::quantile(d$val, c(0.25, 0.5, 0.75), names = FALSE, na.rm = TRUE)
    iqr <- q[3] - q[1]
    lw <- min(d$val[d$val >= q[1] - 1.5 * iqr], na.rm = TRUE)
    uw <- max(d$val[d$val <= q[3] + 1.5 * iqr], na.rm = TRUE)
    data.frame(ymin = lw, lower = q[1], middle = q[2], upper = q[3], ymax = uw)
  }) %>% dplyr::ungroup()

  rng <- bx %>% dplyr::group_by(panel) %>% dplyr::summarise(
    lo0 = min(ymin), hi0 = max(ymax), .groups = "drop") %>%
    dplyr::mutate(span = hi0 - lo0,
                  lo = lo0 - 0.10 * span,     # 아래쪽 여백 (n= 라벨 자리)
                  hi = hi0 + 0.50 * span)     # 위쪽 여백 (검정 결과 주석 자리)
  blank <- rng %>% dplyr::select(panel, lo, hi) %>%
    tidyr::pivot_longer(c(lo, hi), values_to = "val") %>%
    dplyr::mutate(grp = factor("Tumor", levels = levels(E$grp)))

  Ej <- E %>% dplyr::left_join(rng[, c("panel","lo","hi")], by = "panel")
  n_clip <- Ej %>% dplyr::group_by(panel) %>%
    dplyr::summarise(n_outside = sum(val < lo | val > hi), .groups = "drop")
  Ej <- Ej %>% dplyr::filter(val >= lo, val <= hi)
  clip_note <- if (sum(n_clip$n_outside) > 0)
    sprintf("\nVertical axes are set by the whiskers; %d individual points lying outside them are not drawn.\nBox statistics use all samples; the per-panel counts are given in Figure2_expression_source.csv.",
            sum(n_clip$n_outside)) else ""
  message(sprintf("  [Fig2] 표시 범위 밖으로 제외된 점: %d개 (상자 통계량은 전체 자료 사용)",
                  sum(n_clip$n_outside)))

  p2 <- ggplot(bx, aes(x = grp, fill = grp)) +
    geom_blank(data = blank, aes(x = grp, y = val), inherit.aes = FALSE) +
    geom_jitter(data = Ej, aes(x = grp, y = val), inherit.aes = FALSE,
                width = 0.16, size = 0.28, alpha = 0.22, colour = "grey20") +
    geom_boxplot(aes(ymin = ymin, lower = lower, middle = middle,
                     upper = upper, ymax = ymax),
                 stat = "identity", width = 0.55, linewidth = 0.3,
                 colour = "grey25", alpha = 0.9) +
    geom_text(data = nlab, aes(x = grp, y = -Inf, label = paste0("n=", n)),
              inherit.aes = FALSE, vjust = -0.6, size = 2.3, colour = INK["muted"]) +
    geom_text(data = lab, aes(x = 0.55, y = Inf, label = txt), inherit.aes = FALSE,
              hjust = 0, vjust = 1.45, size = 2.3, colour = INK["secondary"], lineheight = 1.05) +
    facet_wrap(~ panel, scales = "free_y", nrow = 2) +
    ## 여러 레이어에서 수준 순서가 뒤집히는 것을 막기 위해 축 순서를 명시한다
    scale_x_discrete(limits = c("GTEx normal", "Adjacent normal", "Tumor")) +
    scale_fill_manual(values = c("GTEx normal" = PAL[["blue"]],
                                 "Adjacent normal" = PAL[["aqua"]],
                                 "Tumor" = PAL[["orange"]])) +
    labs(x = NULL, y = expression(log[2]*"(TPM + 0.001)"),
         title = "Expression regenerated from the raw UCSC Xena TOIL matrix",
         caption = paste0("Boxes show median and interquartile range; whiskers 1.5 x IQR; points are individual samples.\nStatistics are two-sided Wilcoxon rank-sum tests of tumor against each normal reference.",
                          clip_note)) +
    theme_ms() + theme(panel.grid.major.x = element_blank(),
                       panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
                       axis.text.x = element_text(colour = INK["primary"]),
                       legend.position = "none",
                       plot.margin = margin(6, 8, 4, 6))
  save_fig(p2, "Figure2_expression", 180, 150)
  save_src(dplyr::left_join(ann, n_clip, by = "panel"), "Figure2_expression")
}, silent = FALSE)

## =====================================================================
## Figure 3B — 간 vs 폐 공발현 (95% CI 포함)
## =====================================================================
hdr("Figure 3B · 공발현 조직 대비")
try({
  CO  <- rd("12_coexpression_full.csv")
  DIF <- rd("12_tissue_difference.csv")
  ord <- DIF[order(DIF$q_diff), "gene"]
  CO$gene <- factor(CO$gene, levels = rev(ord))
  CO$cohort <- factor(CO$cohort, levels = c("LIHC","LUAD"))
  star <- DIF %>% dplyr::transmute(gene = factor(gene, levels = rev(ord)),
                                   mark = ifelse(q_diff < 0.05, "*", ""))
  p3 <- ggplot(CO, aes(rho, gene, colour = cohort)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey60") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.6,
                   position = position_dodge(width = 0.62)) +
    geom_point(size = 1.9, position = position_dodge(width = 0.62)) +
    geom_text(data = star, aes(x = 0.34, y = gene, label = mark), inherit.aes = FALSE,
              size = 4, colour = INK["primary"], vjust = 0.75) +
    scale_colour_manual(values = c(LIHC = PAL[["orange"]], LUAD = PAL[["blue"]])) +
    ## seq() 는 -0.29999999... 같은 부동소수점 값을 만들어 축 라벨이 깨지므로
    ## 눈금과 라벨을 문자열로 직접 지정한다.
    scale_x_continuous(limits = c(-0.35, 0.38),
                       breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3),
                       labels = c("-0.3","-0.2","-0.1","0","0.1","0.2","0.3")) +
    labs(x = "Spearman correlation with CYB5R3 (95% CI)", y = NULL,
         title = "Exploratory pairwise associations of CYB5R3, liver versus lung",
         caption = "* Benjamini-Hochberg q < 0.05 for the Fisher r-to-z test of the liver-versus-lung difference.\nGenes ordered by that q value. TCGA-LIHC n = 371, TCGA-LUAD n = 515.") +
    theme_ms()
  save_fig(p3, "Figure3B_coexpression", 130, 130)
  save_src(merge(CO, DIF[, c("gene","z_diff","p_diff","q_diff","verdict")], by = "gene"),
           "Figure3B_coexpression")
}, silent = FALSE)

## =====================================================================
## Figure 5 — MASLD 코호트: 코호트별 추정치 + 통합 + I²
## =====================================================================
hdr("Figure 5 · MASLD 조직학")
try({
  TR <- rd("14_masld_trend_all.csv"); MT <- rd("14_masld_meta.csv")
  GEN <- c("CYB5R3","MTARC1","PARP16","SCD","NQO1")
  FEA <- c("fibrosis stage","NAFLD activity score")
  tr <- TR[TR$gene %in% GEN & TR$histology %in% FEA, ]
  mt <- MT[MT$gene %in% GEN & MT$histology %in% FEA, ]
  tr$gene <- factor(tr$gene, levels = rev(GEN)); mt$gene <- factor(mt$gene, levels = rev(GEN))
  tr$histology <- factor(tr$histology, levels = FEA); mt$histology <- factor(mt$histology, levels = FEA)
  mt$i2lab <- sprintf("I² = %.0f%%", mt$I2)

  p5 <- ggplot() +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey60") +
    geom_errorbarh(data = mt, aes(y = gene, xmin = lo, xmax = hi),
                   height = 0, linewidth = 0.7, colour = INK["primary"]) +
    geom_point(data = tr, aes(rho, gene, shape = cohort),
               size = 1.7, colour = PAL[["blue"]], fill = "white", stroke = 0.6,
               position = position_nudge(y = 0.22)) +
    geom_point(data = mt, aes(rho_pooled, gene), size = 2.6,
               shape = 23, fill = PAL[["orange"]], colour = INK["primary"], stroke = 0.4) +
    geom_text(data = mt, aes(x = 0.62, y = gene, label = i2lab),
              hjust = 1, size = 2.3, colour = INK["secondary"]) +
    facet_wrap(~ histology, nrow = 1) +
    scale_shape_manual(values = c(21, 24, 22)) +
    scale_x_continuous(limits = c(-0.55, 0.65),
                       breaks = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6),
                       labels = c("-0.4","-0.2","0","0.2","0.4","0.6")) +
    labs(x = "Spearman correlation with histological grade", y = NULL,
         title = "The two axis members diverge across the human MASLD spectrum",
         subtitle = "Open symbols, individual cohorts; filled diamonds with bars, random-effects pooled estimate and 95% CI",
         caption = "GSE130970 n = 78, GSE135251 n = 216, GSE162694 n = 99-102. I-squared is between-cohort heterogeneity.\nCYB5R3 estimates differ in sign between cohorts; MTARC1, PARP16 and SCD do not.") +
    theme_ms() + theme(legend.position = "top")
  save_fig(p5, "Figure5_masld", 180, 105)
  save_src(merge(tr, mt[, c("gene","histology","rho_pooled","lo","hi","q","I2")],
                 by = c("gene","histology"), all.x = TRUE), "Figure5_masld")
}, silent = FALSE)

## =====================================================================
## Figure 6 — 예후 추정치 forest (TCGA 주분석 + GSE14520 비보정/보정)
## =====================================================================
hdr("Figure 6 · 예후 forest")
try({
  PR <- rd("11_cox_primary.csv")
  GX <- rd("11_gse14520_adjusted.csv")
  pick <- function(df, pat) df[grepl(pat, df$model), ]
  f <- rbind(
    data.frame(cohort = "TCGA-LIHC", adjust = "Adjusted (age, sex, stage)",
               gene = ifelse(pick(PR, "^PRIMARY")$term == "z5", "CYB5R3", "MTARC1"),
               HR = pick(PR, "^PRIMARY")$HR, lo = pick(PR, "^PRIMARY")$CI_low,
               hi = pick(PR, "^PRIMARY")$CI_high, p = pick(PR, "^PRIMARY")$p,
               n = pick(PR, "^PRIMARY")$n, ev = pick(PR, "^PRIMARY")$events),
    data.frame(cohort = "GSE14520", adjust = "Unadjusted",
               gene = ifelse(grepl("CYB5R3", pick(GX, "UNADJUSTED")$model), "CYB5R3", "MTARC1"),
               HR = pick(GX, "UNADJUSTED")$HR, lo = pick(GX, "UNADJUSTED")$CI_low,
               hi = pick(GX, "UNADJUSTED")$CI_high, p = pick(GX, "UNADJUSTED")$p,
               n = pick(GX, "UNADJUSTED")$n, ev = pick(GX, "UNADJUSTED")$events),
    data.frame(cohort = "GSE14520", adjust = "Adjusted (age, sex, TNM)",
               gene = ifelse(grepl("CYB5R3", pick(GX, "ADJUSTED \\(age")$model), "CYB5R3", "MTARC1"),
               HR = pick(GX, "ADJUSTED \\(age")$HR, lo = pick(GX, "ADJUSTED \\(age")$CI_low,
               hi = pick(GX, "ADJUSTED \\(age")$CI_high, p = pick(GX, "ADJUSTED \\(age")$p,
               n = pick(GX, "ADJUSTED \\(age")$n, ev = pick(GX, "ADJUSTED \\(age")$events))
  f$row <- factor(paste(f$cohort, f$adjust, sep = "\n"),
                  levels = rev(unique(paste(f$cohort, f$adjust, sep = "\n"))))
  f$lab <- sprintf("%.2f (%.2f-%.2f), p = %s", f$HR, f$lo, f$hi, fmt_p(f$p))

  p6 <- ggplot(f, aes(HR, row, colour = gene)) +
    geom_vline(xintercept = 1, linewidth = 0.35, colour = "grey60") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.7,
                   position = position_dodge(width = 0.55)) +
    geom_point(size = 2.2, position = position_dodge(width = 0.55)) +
    ## 수치 라벨은 자료 영역 오른쪽의 전용 구역에 배치한다
    ## (이전 판에서는 신뢰구간 끝과 라벨이 겹쳤다).
    geom_text(aes(x = 2.55, label = lab), hjust = 0, size = 2.3,
              position = position_dodge(width = 0.55), show.legend = FALSE) +
    scale_colour_manual(values = c(CYB5R3 = PAL[["orange"]], MTARC1 = PAL[["blue"]])) +
    scale_x_log10(limits = c(0.5, 5.8), breaks = c(0.6, 0.8, 1.0, 1.25, 1.6, 2.0),
                  labels = c("0.60","0.80","1.00","1.25","1.60","2.00"),
                  expand = expansion(mult = c(0.01, 0))) +
    labs(x = "Hazard ratio for overall survival per SD of expression (log scale)", y = NULL,
         title = "Prognostic estimates in the primary and external cohorts",
         caption = "TCGA-LIHC n = 339, 114 deaths. GSE14520 n = 242, 96 deaths (unadjusted) and n = 225, 86 deaths (adjusted).\nAdjustment attenuates both external estimates below conventional significance.") +
    theme_ms() + theme(axis.text.y = element_text(colour = INK["primary"], size = 8))
  save_fig(p6, "Figure6_prognostic", 180, 85)
  save_src(f, "Figure6_prognostic")
}, silent = FALSE)

## =====================================================================
## Supplementary Figure S6 — 환자 흐름
## =====================================================================
hdr("Figure S6 · 환자 흐름")
try({
  FL <- rd("11_flow.csv")
  FL <- FL[!grepl("Median follow-up", FL$step), ]
  FL$order <- seq_len(nrow(FL))
  ## 철자 표기는 미국식으로 통일 (CSV 가 이전 판이어도 그림은 올바르게 나오도록)
  FL$step <- gsub("tumours", "tumors", gsub("Tumour", "Tumor", FL$step))
  FL$label <- sprintf("%s\nn = %s", FL$step, format(FL$n, big.mark = ""))
  ps6 <- ggplot(FL, aes(x = 1, y = -order)) +
    geom_tile(width = 0.9, height = 0.62, fill = "grey96", colour = "grey70", linewidth = 0.3) +
    geom_text(aes(label = label), size = 2.6, colour = INK["primary"], lineheight = 1.1) +
    geom_segment(data = FL[-nrow(FL), ],
                 aes(x = 1, xend = 1, y = -order - 0.31, yend = -order - 0.69),
                 arrow = arrow(length = unit(1.6, "mm"), type = "closed"),
                 linewidth = 0.35, colour = "grey55") +
    scale_x_continuous(limits = c(0.4, 1.6)) +
    labs(title = "Patient flow, TCGA-LIHC primary analysis",
         caption = "Each step reports the number of patients remaining. The primary analysis set is the final row.") +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 10, colour = INK["primary"]),
          plot.caption = element_text(size = 7, colour = INK["muted"], hjust = 0),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(8, 8, 8, 8))
  save_fig(ps6, "FigureS6_patient_flow", 120, 130)
}, silent = FALSE)

## =====================================================================
## Supplementary Figure S7 — Schoenfeld 잔차 (주모형 재적합)
## =====================================================================
hdr("Figure S7 · Schoenfeld 잔차")
try({
  ## 11_survival_primary.R 의 주분석 집합(n = 339)을 그대로 재구성한다.
  ## (CYB5R3 → 원발종양 1인1검체 → MTARC1 병합 → 생존 → 완전한 공변량)
  c5 <- fetch_gene(XH$tcga, XDS$lihc_hiseq, "CYB5R3"); stopifnot(!is.null(c5))
  v5 <- tovec(c5); names(v5) <- nb(names(v5))
  C <- data.frame(sample = names(v5), CYB5R3 = as.numeric(v5))
  C <- C[typ(C$sample) %in% TUMOR_CODES, ]; C$patient <- pat(C$sample)
  C <- C %>% dplyr::group_by(patient) %>% dplyr::slice(1) %>% dplyr::ungroup()

  mt <- fetch_gene(XH$toil, XDS$toil_tpm, "MTARC1"); stopifnot(!is.null(mt))
  vm <- tovec(mt); names(vm) <- nb(names(vm))
  M <- data.frame(sample = names(vm), MTARC1 = as.numeric(vm))
  M <- M[grepl("^TCGA", M$sample) & typ(M$sample) %in% TUMOR_CODES, ]
  M$patient <- pat(M$sample)
  M <- M %>% dplyr::group_by(patient) %>% dplyr::slice(1) %>% dplyr::ungroup()
  E <- dplyr::inner_join(C[, c("sample","patient","CYB5R3")],
                         M[, c("patient","MTARC1")], by = "patient")

  xq <- UCSCXenaTools::XenaGenerate() %>% UCSCXenaTools::XenaFilter(filterDatasets = XDS$lihc_clin)
  xd <- UCSCXenaTools::XenaQuery(xq) %>%
        UCSCXenaTools::XenaDownload(destdir = CACHE, trans_slash = TRUE, force = FALSE)
  cl <- UCSCXenaTools::XenaPrepare(xd); if (is.list(cl) && !is.data.frame(cl)) cl <- cl[[1]]
  cl <- as.data.frame(cl); names(cl)[1] <- "sample"; cl$patient <- pat(cl$sample)
  gk <- function(p) grep(p, names(cl), ignore.case = TRUE, value = TRUE)[1]
  cl$age <- suppressWarnings(as.numeric(as.character(cl[[gk("^age_at_initial|^age$")]])))
  cl$sex <- as.character(cl[[gk("^gender$|^sex$")]])
  sr <- as.character(cl[[gk("pathologic_stage|^stage$|ajcc.*stage")]])
  ## 11번 스크립트와 동일한 병기 부호화
  cl$stage_bin <- dplyr::case_when(
    stringr::str_detect(sr, "Stage III|Stage IV")                  ~ "III-IV",
    stringr::str_detect(sr, "Stage I{1,2}$|Stage I |Stage II")      ~ "I-II",
    stringr::str_detect(sr, "^Stage I$")                            ~ "I-II",
    TRUE ~ NA_character_)
  cl <- cl %>% dplyr::group_by(patient) %>% dplyr::slice(1) %>% dplyr::ungroup()

  sv <- UCSCXenaTools::fetch_dense_values(XH$pan, XDS$pan_surv, c("OS","OS.time"), use_probeMap = FALSE)
  sv <- as.data.frame(t(sv)); sv$patient <- pat(rownames(sv))
  sv$OS <- suppressWarnings(as.numeric(as.character(sv$OS)))
  sv$OS.time <- suppressWarnings(as.numeric(as.character(sv$OS.time)))
  sv <- sv %>% dplyr::group_by(patient) %>% dplyr::slice(1) %>% dplyr::ungroup()

  d <- E %>% dplyr::left_join(cl[, c("patient","age","sex","stage_bin")], by = "patient") %>%
             dplyr::left_join(sv[, c("patient","OS","OS.time")], by = "patient") %>%
             dplyr::filter(!is.na(OS), !is.na(OS.time), OS.time > 0) %>%
             dplyr::filter(!is.na(age), !is.na(sex), !is.na(stage_bin))
  d$z5 <- as.numeric(scale(d$CYB5R3)); d$sex <- factor(d$sex)
  d$stage_bin <- factor(d$stage_bin, c("I-II","III-IV"))
  message(sprintf("  [S7] 주분석 집합 재구성: n = %d (사망 %d) — 11번 스크립트와 일치해야 함(339 / 114)",
                  nrow(d), sum(d$OS == 1)))
  if (nrow(d) != 339)
    warning(sprintf("Figure S7 의 분석집합 n = %d 이며 주분석(339)과 다릅니다. 본문 수치와 대조하십시오.", nrow(d)))
  fit <- survival::coxph(survival::Surv(OS.time, OS) ~ z5 + age + sex + stage_bin, data = d)
  zph <- survival::cox.zph(fit)
  ## Schoenfeld 패널 제목을 사람이 읽는 이름으로
  s7lab <- c(z5 = "CYB5R3 (per SD)", age = "Age (years)", sex = "Sex",
             stage_bin = "AJCC stage III-IV vs I-II")
  s7nm <- function(k) if (k %in% names(s7lab)) unname(s7lab[k]) else k
  gp <- zph$table["GLOBAL", "p"]
  message(sprintf("  [S7] 전역 Schoenfeld p = %.3f", gp))
  w_res(data.frame(term = rownames(zph$table), chisq = zph$table[, "chisq"],
                   df = zph$table[, "df"], p = zph$table[, "p"],
                   n = nrow(d), deaths = sum(d$OS == 1)),
        "17_figS7_ph_test.csv")
  png(file.path(FIG, "FigureS7_schoenfeld.png"), width = MM(180)*300, height = MM(150)*300,
      res = 300, bg = "white")
  par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), cex = 0.8, col.axis = "grey35", col.lab = "black")
  for (i in seq_len(nrow(zph$table) - 1)) {
    plot(zph[i], col = PAL[["orange"]], lwd = 1.6, resid = TRUE,
         main = sprintf("%s   (p = %s)", s7nm(rownames(zph$table)[i]), fmt_p(zph$table[i, "p"])),
         ylab = sprintf("Beta(t) for %s", s7nm(rownames(zph$table)[i])))
    abline(h = 0, col = "grey60", lty = 2)
  }
  dev.off()
  tiff(file.path(FIG, "FigureS7_schoenfeld.tiff"), width = MM(180)*300, height = MM(150)*300,
       res = 300, compression = "lzw", bg = "white")
  par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), cex = 0.8)
  for (i in seq_len(nrow(zph$table) - 1)) {
    plot(zph[i], col = PAL[["orange"]], lwd = 1.6, resid = TRUE,
         main = sprintf("%s   (p = %s)", s7nm(rownames(zph$table)[i]), fmt_p(zph$table[i, "p"])),
         ylab = sprintf("Beta(t) for %s", s7nm(rownames(zph$table)[i])))
    abline(h = 0, col = "grey60", lty = 2)
  }
  dev.off()
  message("  → FigureS7_schoenfeld.png / .tiff  (global p = ", fmt_p(zph$table["GLOBAL","p"]), ")")
}, silent = FALSE)

save_session("17_figures")
message("\n[17] 완료. figures/ 폴더를 확인하십시오.")
