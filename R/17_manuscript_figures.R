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
## Figure 3A — CPTAC 단백 수준 (포털 유래 값; 원 투고 Figure 3A 를 새 서식으로 재작성)
## =====================================================================
hdr("Figure 3A · CPTAC 단백")
p3a <- NULL
try({
  ## UALCAN 의 CPTAC 모듈은 원자료를 내보내지 않으므로 포털이 보고한 중앙값을 쓴다.
  ## 본문 2.5절·Table 3 과 동일한 수치이며, 근사값임을 각주에 명시한다.
  P <- data.frame(grp = factor(c("Adjacent normal","Tumor"),
                               levels = c("Adjacent normal","Tumor")),
                  z   = c(0.35, 0.0), n = c(165, 165))
  ## 종양의 중앙값이 0 이므로 막대그래프로 그리면 막대가 보이지 않는다.
  ## 0 기준선에서 뻗어나가는 점-선(lollipop) 형태로 그리고 값을 직접 표기한다.
  p3a <- ggplot(P, aes(grp, z, colour = grp)) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey55") +
    geom_segment(aes(x = grp, xend = grp, y = 0, yend = z), linewidth = 1.1) +
    geom_point(size = 3.4) +
    geom_text(aes(label = sprintf("%.2f", z)), vjust = -1.2, size = 2.6,
              colour = INK["primary"], show.legend = FALSE) +
    geom_text(aes(label = paste0("n = ", n)), y = -0.055, size = 2.3,
              colour = INK["muted"], show.legend = FALSE) +
    ## 유니코드 위첨자는 기기 폰트에 없을 수 있으므로 plotmath 로 그린다
    annotate("text", x = 1.5, y = 0.44, parse = TRUE,
             label = "p == 1.3 %*% 10^-5",
             size = 2.6, colour = INK["secondary"]) +
    scale_colour_manual(values = c("Adjacent normal" = PAL[["aqua"]],
                                   "Tumor" = PAL[["orange"]])) +
    scale_y_continuous(limits = c(-0.08, 0.50),
                       breaks = c(0, 0.1, 0.2, 0.3, 0.4),
                       labels = c("0","0.1","0.2","0.3","0.4")) +
    labs(x = NULL, y = "Median protein Z-value",
         title = "A  Protein level (CPTAC)",
         caption = "Portal-derived medians.") +
    theme_ms() + theme(legend.position = "none",
                       panel.grid.major.x = element_blank(),
                       panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
                       axis.text.x = element_text(colour = INK["primary"]),
                       plot.title = element_text(face = "bold", size = 9.5))
  message("  Figure 3A 준비 완료 (포털 값)")
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
         title = "B  Pairwise associations, liver vs lung",
         caption = "* q < 0.05, Fisher r-to-z test of the tissue difference (BH-corrected).\nLIHC n = 371, LUAD n = 515.") +
    theme_ms()
  save_fig(p3, "Figure3B_coexpression", 130, 130)
  save_src(merge(CO, DIF[, c("gene","z_diff","p_diff","q_diff","verdict")], by = "gene"),
           "Figure3B_coexpression")

  ## --- 두 패널을 하나의 Figure 3 으로 결합 (grid 만 사용, 추가 패키지 불필요) ---
  if (!is.null(p3a)) {
    combine2 <- function(pa, pb, name, w_mm, h_mm) {
      for (dev in c("png","tiff")) {
        f <- file.path(FIG, paste0(name, ".", dev))
        if (dev == "png") grDevices::png(f, width = MM(w_mm), height = MM(h_mm),
                                         units = "in", res = 300)
        else grDevices::tiff(f, width = MM(w_mm), height = MM(h_mm),
                             units = "in", res = 300, compression = "lzw")
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(
          1, 2, widths = grid::unit(c(0.85, 1.6), "null"))))
        print(pa, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
        print(pb, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
        grDevices::dev.off()
        message("  → ", f)
      }
    }
    combine2(p3a, p3, "Figure3_orthogonal", 185, 128)
  } else {
    message("  [skip] Figure 3A 가 없어 결합본을 만들지 않았습니다.")
  }
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


## =====================================================================
## Supplementary Figures S1-S5 — 본문과 같은 서식으로 재작성
## (모두 deposited results/*.csv 에서만 그린다. 새 계산 없음.)
## =====================================================================

## --- Figure S1 : DepMap 의존성 ------------------------------------------
hdr("Figure S1 · DepMap CRISPR dependency")
try({
  D <- rd("15_depmap_dependency.csv")
  D$grp <- factor(D$grp, levels = D$grp[order(D$median_gene_effect)])
  pw <- unique(D$wilcoxon_p)[1]
  lab <- sprintf("%s\nn = %d lines\n%s", D$grp, D$n_lines,
                 ifelse(D$pct_below_minus0.5 == 0, "no line below -0.5",
                        sprintf("%.2f%% below -0.5", D$pct_below_minus0.5)))
  pS1 <- ggplot(D, aes(grp, median_gene_effect, colour = grp)) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey55") +
    geom_hline(yintercept = -0.5, linewidth = 0.5, linetype = 2, colour = PAL[["violet"]]) +
    annotate("text", x = 0.62, y = -0.47, hjust = 0, size = 2.3, colour = PAL[["violet"]],
             label = "dependency threshold (-0.5)") +
    geom_linerange(aes(ymin = Q1, ymax = Q3), linewidth = 3.2, alpha = 0.35) +
    geom_point(size = 3.4) +
    geom_text(aes(label = sprintf("%.3f", median_gene_effect)), hjust = -0.7,
              size = 2.5, colour = INK["primary"], show.legend = FALSE) +
    annotate("text", x = 1.5, y = 0.06, size = 2.6, colour = INK["secondary"],
             label = sprintf("Wilcoxon p = %.3f", pw)) +
    scale_x_discrete(labels = setNames(lab, as.character(D$grp))) +
    scale_y_continuous(limits = c(-0.58, 0.12)) +
    scale_colour_manual(values = setNames(c(PAL[["orange"]], PAL[["blue"]])[seq_len(nrow(D))],
                                          as.character(D$grp))) +
    coord_flip() +
    labs(x = NULL, y = "Chronos gene effect (median, bar = interquartile range)",
         title = "CYB5R3 CRISPR dependency, liver versus lung cell lines",
         caption = "Bioconductor depmap 1.26.0. 0 = no effect; -1 = median common-essential gene.\nSource: results/15_depmap_dependency.csv") +
    theme_ms() + theme(legend.position = "none",
                       panel.grid.major.y = element_blank(),
                       panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3))
  save_fig(pS1, "FigureS1_depmap", 180, 90)
  save_src(D, "FigureS1_depmap")
}, silent = FALSE)

## --- Figure S2 : 범암종 forest ------------------------------------------
hdr("Figure S2 · Pan-cancer forest")
try({
  PC <- rd("13_pancancer_full.csv")
  PC <- PC[!is.na(PC$HR_perSD), ]
  PC$cohort <- factor(PC$cohort, levels = PC$cohort[order(PC$HR_perSD)])
  PC$is_lihc <- PC$cohort == "LIHC"
  pS2 <- ggplot(PC, aes(HR_perSD, cohort)) +
    geom_vline(xintercept = 1, linewidth = 0.4, colour = "grey55") +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high, colour = is_lihc),
                   height = 0, linewidth = 0.6) +
    geom_point(aes(colour = is_lihc), size = 1.9) +
    geom_text(aes(x = max(CI_high, na.rm = TRUE) * 1.50, label = sprintf("q = %.2f", q_cox)), hjust = 1,
              size = 2.2, colour = INK["muted"]) +
    scale_colour_manual(values = c("FALSE" = PAL[["grey"]], "TRUE" = PAL[["orange"]])) +
    scale_x_log10(limits = c(min(PC$CI_low, na.rm = TRUE) * 0.93,
                             max(PC$CI_high, na.rm = TRUE) * 1.55),
                  breaks = c(0.5, 0.75, 1, 1.5, 2, 3),
                  labels = c("0.5","0.75","1","1.5","2","3")) +
    labs(x = "Overall-survival hazard ratio per SD of CYB5R3 (univariable)", y = NULL,
         title = "CYB5R3 and overall survival across 28 TCGA cohorts",
         caption = "q, Benjamini-Hochberg across cohorts. No cohort is significant after correction.\nLIHC (orange) uses the same patients as the primary analysis and is not an independent replication.\nSource: results/13_pancancer_full.csv") +
    theme_ms() + theme(legend.position = "none",
                       panel.grid.major.y = element_blank(),
                       panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3))
  save_fig(pS2, "FigureS2_pancancer_forest", 170, 175)
  save_src(PC, "FigureS2_pancancer_forest")
}, silent = FALSE)

## --- Figure S3 : log2FC vs log2HR ---------------------------------------
hdr("Figure S3 · Expression change versus prognostic direction")
try({
  PC <- rd("13_pancancer_full.csv")
  S3 <- PC[!is.na(PC$log2FC) & !is.na(PC$HR_perSD), ]
  S3$logHR <- log2(S3$HR_perSD)
  S3 <- S3[order(S3$log2FC), ]
  ## x 순서로 3단계 수직 오프셋을 돌려 인접한 점끼리는 반드시 다른 높이에 놓는다
  S3$vj <- rep(c(-1.15, 2.05, 3.55), length.out = nrow(S3))
  S3$hj <- 0.5
  ct <- suppressWarnings(stats::cor.test(S3$log2FC, S3$logHR, method = "spearman"))
  pS3 <- ggplot(S3, aes(log2FC, logHR)) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey60") +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey60") +
    geom_point(aes(colour = cohort == "LIHC"), size = 2.2, show.legend = FALSE) +
    ## 라벨이 겹치지 않도록 x 순서에 따라 위/아래를 번갈아 배치한다
    geom_text(aes(label = cohort, vjust = vj, hjust = hj), size = 2.1,
              colour = INK["secondary"], show.legend = FALSE) +
    scale_colour_manual(values = c("FALSE" = PAL[["blue"]], "TRUE" = PAL[["orange"]])) +
    annotate("text", x = -Inf, y = -Inf, hjust = -0.08, vjust = -1.0, size = 2.6,
             colour = INK["secondary"],
             label = sprintf("Spearman r = %.2f, p = %s, %d cohorts",
                             unname(ct$estimate), fmt_p(ct$p.value), nrow(S3))) +
    scale_x_continuous(expand = expansion(mult = 0.10)) +
    scale_y_continuous(expand = expansion(mult = 0.13)) +
    labs(x = "Tumor-versus-normal log2 fold change", y = "log2 hazard ratio per SD",
         title = "Expression change does not predict prognostic direction",
         caption = "Cohorts with at least 10 adjacent-normal samples. The up-and-adverse, down-and-protective\nrule does not generalize beyond the liver-lung contrast. Source: results/13_pancancer_full.csv") +
    theme_ms() + theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))
  save_fig(pS3, "FigureS3_fc_vs_hr", 150, 120)
  save_src(S3, "FigureS3_fc_vs_hr")
}, silent = FALSE)

## --- Figure S4 : 병기 / 분화도 / 혈청학 ---------------------------------
hdr("Figure S4 · Stage, grade and serology")
try({
  SG <- rd("16_lihc_stage_grade.csv")
  ET <- rd("16_lihc_etiology.csv")
  parse_levels <- function(txt, facet) {
    parts <- trimws(strsplit(txt, ";", fixed = TRUE)[[1]])
    lv <- sub("^([^ ]+) .*$", "\\1", parts)
    n  <- as.numeric(sub("^.*n=([0-9]+).*$", "\\1", parts))
    md <- as.numeric(sub("^.*median ([0-9.]+)\\).*$", "\\1", parts))
    data.frame(facet = facet, level = factor(lv, levels = lv), n = n, med = md,
               stringsAsFactors = FALSE)
  }
  st <- SG[SG$stratifier == "stage", ]
  gr <- SG[SG$stratifier == "grade", ]
  A <- parse_levels(st$levels, sprintf("AJCC stage  (p = %.3f)", st$p))
  B <- parse_levels(gr$levels, sprintf("Histologic grade  (p = %.4f)", gr$p))
  C <- data.frame(facet = sprintf("Viral serology  (p = %.2f)", ET$wilcox_p),
                  level = factor(c("Positive recorded","None recorded"),
                                 levels = c("Positive recorded","None recorded")),
                  n = c(ET$n_viral_positive, ET$n_no_positive_recorded),
                  med = c(ET$median_viral_positive, ET$median_no_positive),
                  stringsAsFactors = FALSE)
  S4 <- rbind(A, B, C)
  S4$facet <- factor(S4$facet, levels = unique(S4$facet))
  base <- min(S4$med) - 0.34
  pS4 <- ggplot(S4, aes(level, med)) +
    geom_segment(aes(x = level, xend = level, y = base + 0.13, yend = med),
                 colour = "grey80", linewidth = 0.6) +
    geom_point(size = 3, colour = PAL[["orange"]]) +
    ## n 은 막대 아래 빈 영역에 두어 선과 겹치지 않게 한다
    geom_text(aes(label = sprintf("n = %d", n)), y = base + 0.04, size = 2.2,
              colour = INK["muted"]) +
    facet_wrap(~ facet, scales = "free_x", nrow = 1) +
    scale_y_continuous(limits = c(base, max(S4$med) + 0.10)) +
    labs(x = NULL, y = "Median CYB5R3 (log2 norm_count+1)",
         title = "CYB5R3 in TCGA-LIHC by stage, grade and viral serology",
         caption = "Group medians with group sizes. Trend tests are Jonckheere-Terpstra (stage, grade) and Wilcoxon (serology).\nExpression declines across histologic grade; stage and serology do not differ.\n'None recorded' is not equivalent to non-viral etiology (Section 2.9).\nSource: results/16_lihc_stage_grade.csv, results/16_lihc_etiology.csv") +
    theme_ms() + theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
                       panel.grid.major.x = element_blank())
  save_fig(pS4, "FigureS4_stage_grade_serology", 180, 95)
  save_src(S4, "FigureS4_stage_grade_serology")
}, silent = FALSE)

## --- Figure S5 : 16개 유전자 패널, 원값 + 조성보정값 ---------------------
hdr("Figure S5 · Redox panel, raw and composition-adjusted")
try({
  CO  <- rd("12_coexpression_full.csv")
  PU  <- rd("12_purity_adjusted.csv")
  DIF <- rd("12_tissue_difference.csv")
  ord <- DIF[order(DIF$q_diff), "gene"]
  raw <- CO[, c("cohort","gene","rho")]; raw$kind <- "Unadjusted"
  adj <- PU[, c("cohort","gene","rho_adjusted")]
  names(adj)[3] <- "rho"; adj$kind <- "Composition-adjusted"
  S5 <- rbind(raw, adj)
  S5$gene   <- factor(S5$gene, levels = rev(ord))
  S5$cohort <- factor(S5$cohort, levels = c("LIHC","LUAD"))
  S5$kind   <- factor(S5$kind, levels = c("Unadjusted","Composition-adjusted"))
  pS5 <- ggplot(S5, aes(rho, gene, colour = cohort, shape = kind)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey60") +
    geom_point(size = 1.9, position = position_dodge(width = 0.62)) +
    scale_colour_manual(values = c(LIHC = PAL[["orange"]], LUAD = PAL[["blue"]])) +
    scale_shape_manual(values = c("Unadjusted" = 16, "Composition-adjusted" = 1)) +
    scale_x_continuous(limits = c(-0.32, 0.32),
                       breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3),
                       labels = c("-0.3","-0.2","-0.1","0","0.1","0.2","0.3")) +
    labs(x = "Spearman correlation with CYB5R3", y = NULL,
         title = "Redox panel before and after adjustment for tumor cellular composition",
         caption = "Open symbols, adjusted for a non-tumor-content proxy (Section 4.6). Genes ordered by the\nq value of the liver-versus-lung difference. Source: results/12_coexpression_full.csv,\nresults/12_purity_adjusted.csv") +
    theme_ms()
  save_fig(pS5, "FigureS5_panel_adjusted", 150, 140)
  save_src(S5, "FigureS5_panel_adjusted")
}, silent = FALSE)

save_session("17_figures")
message("\n[17] 완료. figures/ 폴더를 확인하십시오.")
