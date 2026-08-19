# =====================================================================
#  12_coexpression_tissue.R
#  목적 (심사 대응):
#   · R1-3 / R3-2 "Formal comparisons of correlation coefficients between LIHC
#     and LUAD" → Fisher r-to-z 검정 + 차이의 95% CI
#   · R3 "no confidence intervals … does not adjust for tumour purity …
#     never carries out a statistical test of whether the liver and lung
#     correlations genuinely differ" → 상관계수 95% CI, 종양순도 보정 편상관,
#     순열(permutation) 귀무분포
#   · R3 "Report the full 16-gene co-expression panel with q-values in the main
#     text, including the CAT and SOD2 reversals"
#   · R1 minor 3 "co-expression module 은 과한 표현" → 결과를 exploratory
#     pairwise association 으로 기술 (표 헤더에 반영)
#
#  산출: results/12_coexpression_full.csv        전체 패널 · r · 95%CI · p · q
#        results/12_tissue_difference.csv        Fisher r-to-z 검정 (LIHC vs LUAD)
#        results/12_purity_adjusted.csv          종양순도 보정 편상관
#        results/12_permutation_reversals.csv    부호역전 개수의 순열 귀무분포
#        figures/12_coexpression_liver_lung.png
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
hdr("12 · CYB5R3 공발현 — 조직 대비, 순도 보정, 형식적 검정")

PANEL <- unique(c(REDOX_PANEL, "MTARC1"))
N_PERM <- 10000

## =====================================================================
## 1. 발현 행렬 (원발종양만)
## =====================================================================
get_cohort <- function(co) {
  ds <- XDS[[paste0(tolower(co), "_hiseq")]]
  df <- fetch_matrix(XH$tcga, ds, unique(c("CYB5R3", PANEL, STROMAL_IMMUNE)))
  if (is.null(df)) stop("발현 조회 실패: ", co)
  df$sample <- nb(df$sample)
  df <- df[typ(df$sample) %in% TUMOR_CODES, ]
  df$patient <- pat(df$sample)
  df <- df %>% dplyr::group_by(patient) %>% dplyr::slice(1) %>% dplyr::ungroup()
  message(sprintf("  %s: n=%d primary tumors, genes retrieved = %d",
                  co, nrow(df), sum(c("CYB5R3", PANEL) %in% names(df))))
  as.data.frame(df)
}
EX <- list(LIHC = get_cohort("LIHC"), LUAD = get_cohort("LUAD"))

## =====================================================================
## 2. 종양순도 — TCGA ABSOLUTE (Xena PanCanAtlas). 실패 시 ESTIMATE 대체.
## =====================================================================
hdr("2. 종양순도 추정치 확보")
purity <- NULL
try({
  pp <- xena_table("mastercalls.*abs_tables", exact = XDS$pan_purity, hub = "pancanAtlasHub")
  if (!is.null(pp)) {
    message("  [diag] purity 표 열이름: ", paste(utils::head(names(pp), 15), collapse = " | "))
    scol <- grep("^array$|^sample$|barcode|^aliquot", names(pp), ignore.case = TRUE, value = TRUE)[1]
    pcol <- grep("^purity$|^ABSOLUTE.*purity|^CPE$|^cpe$", names(pp), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(scol) && !is.na(pcol)) {
      purity <- data.frame(sample15 = bc15(pp[[scol]]),
                           purity = suppressWarnings(as.numeric(pp[[pcol]])))
      purity <- purity[!is.na(purity$purity), ]
      purity <- purity[!duplicated(purity$sample15), ]
      message("  ABSOLUTE purity: n = ", nrow(purity))
    }
  }
}, silent = TRUE)
if (is.null(purity) || !nrow(purity))
  message("  ABSOLUTE 순도를 확보하지 못했습니다 → 기질/면역 표지자 기반 대리지표를 사용합니다.\n",
          "  (원고 Methods 에 '대리지표 사용'으로 정직하게 기술하면 심사상 문제 없습니다.)")
## 마지막 대체: 비종양(기질/면역) 함량 대리지표
## — 잘 확립된 기질·면역 표지자의 첫 주성분(부호 반전). 순도의 역방향 대리변수.
proxy_purity <- function(df) {
  gs <- intersect(STROMAL_IMMUNE, names(df))
  if (length(gs) < 5) return(rep(NA_real_, nrow(df)))
  m <- scale(as.matrix(df[, gs, drop = FALSE]))
  m[is.na(m)] <- 0
  -prcomp(m)$x[, 1]     # 부호를 '순도 높음 = 큰 값' 방향으로
}

attach_purity <- function(df) {
  df$sample15 <- bc15(df$sample)
  if (!is.null(purity) && nrow(purity))
    df <- dplyr::left_join(df, purity, by = "sample15") else df$purity <- NA_real_
  df$purity_proxy <- proxy_purity(df)
  df
}
EX <- lapply(EX, attach_purity)
for (co in names(EX))
  message(sprintf("  %s: ABSOLUTE purity available for %d / %d tumors",
                  co, sum(!is.na(EX[[co]]$purity)), nrow(EX[[co]])))

## =====================================================================
## 3. 전체 패널 상관 (r, 95% CI, p, BH q) — 코호트별
## =====================================================================
hdr("3. 전체 패널 Spearman 상관 + 95% CI + BH q")
co_tab <- list()
for (co in names(EX)) {
  d <- EX[[co]]
  gs <- intersect(PANEL, names(d))
  rows <- lapply(gs, function(g) {
    r <- spearman_ci(d$CYB5R3, d[[g]])
    cbind(cohort = co, gene = g, r)
  })
  tb <- dplyr::bind_rows(rows)
  tb$q <- p.adjust(tb$p, method = "BH")     # 패널 내 BH 보정
  co_tab[[co]] <- tb
}
CO <- dplyr::bind_rows(co_tab)
CO$significant_after_FDR <- CO$q < 0.05
print(CO[order(CO$cohort, CO$q), c("cohort","gene","n","rho","lo","hi","p","q")], digits = 3)
w_res(CO, "12_coexpression_full.csv")

## =====================================================================
## 4. 간 vs 폐 상관계수 차이의 형식적 검정 (Fisher r-to-z)
## =====================================================================
hdr("4. LIHC vs LUAD — Fisher r-to-z 검정")
L <- co_tab$LIHC; U <- co_tab$LUAD
mg <- merge(L[, c("gene","n","rho","lo","hi","p","q")],
            U[, c("gene","n","rho","lo","hi","p","q")],
            by = "gene", suffixes = c("_LIHC","_LUAD"))
fz <- fisher_z_diff(mg$rho_LIHC, mg$n_LIHC, mg$rho_LUAD, mg$n_LUAD)
DIF <- cbind(mg, fz)
DIF$q_diff <- p.adjust(DIF$p_diff, method = "BH")
DIF$sign_reversal <- sign(DIF$rho_LIHC) != sign(DIF$rho_LUAD)
DIF$verdict <- dplyr::case_when(
  DIF$q_diff < 0.05 &  DIF$sign_reversal ~ "significant reversal",
  DIF$q_diff < 0.05 & !DIF$sign_reversal ~ "significant difference in magnitude only",
  DIF$sign_reversal                      ~ "apparent reversal, NOT significant after FDR",
  TRUE                                   ~ "no significant tissue difference")
print(DIF[order(DIF$q_diff), c("gene","rho_LIHC","q_LIHC","rho_LUAD","q_LUAD",
                               "z_diff","p_diff","q_diff","verdict")], digits = 3)
w_res(DIF, "12_tissue_difference.csv")

## =====================================================================
## 5. 종양순도 보정 편상관
## =====================================================================
hdr("5. 종양순도 보정 편(partial) Spearman")
pur_rows <- list()
for (co in names(EX)) {
  d <- EX[[co]]
  use_abs <- sum(!is.na(d$purity)) >= 100
  zvar <- if (use_abs) "purity" else "purity_proxy"
  gs <- intersect(PANEL, names(d))
  for (g in gs) {
    raw <- spearman_ci(d$CYB5R3, d[[g]])
    adj <- partial_spearman(d$CYB5R3, d[[g]], d[, zvar, drop = FALSE])
    pur_rows[[length(pur_rows)+1]] <- data.frame(
      cohort = co, gene = g, purity_variable = zvar,
      n_raw = raw$n, rho_raw = raw$rho, p_raw = raw$p,
      n_adj = adj$n, rho_adjusted = adj$rho, p_adjusted = adj$p,
      attenuation = raw$rho - adj$rho)
  }
}
PUR <- dplyr::bind_rows(pur_rows)
PUR <- PUR %>% dplyr::group_by(cohort) %>%
  dplyr::mutate(q_adjusted = p.adjust(p_adjusted, method = "BH")) %>% dplyr::ungroup()
print(as.data.frame(PUR)[, c("cohort","gene","purity_variable","rho_raw","rho_adjusted",
                             "p_adjusted","q_adjusted")], digits = 3)
w_res(PUR, "12_purity_adjusted.csv")

## 순도 보정 후에도 조직 역전이 남는가?
PL <- PUR[PUR$cohort == "LIHC", ]; PU <- PUR[PUR$cohort == "LUAD", ]
mg2 <- merge(PL[, c("gene","n_adj","rho_adjusted")], PU[, c("gene","n_adj","rho_adjusted")],
             by = "gene", suffixes = c("_LIHC","_LUAD"))
fz2 <- fisher_z_diff(mg2$rho_adjusted_LIHC, mg2$n_adj_LIHC,
                     mg2$rho_adjusted_LUAD, mg2$n_adj_LUAD)
DIF2 <- cbind(mg2, fz2); DIF2$q_diff <- p.adjust(DIF2$p_diff, method = "BH")
DIF2$sign_reversal <- sign(DIF2$rho_adjusted_LIHC) != sign(DIF2$rho_adjusted_LUAD)
w_res(DIF2, "12_tissue_difference_purity_adjusted.csv")
print(DIF2[order(DIF2$q_diff), ], digits = 3)

## =====================================================================
## 6. 순열 귀무분포 — 부호역전 개수가 우연 이상인가?
## =====================================================================
hdr("6. 부호역전 개수의 순열 검정")
obs_rev <- sum(DIF$sign_reversal, na.rm = TRUE)
perm_rev <- integer(N_PERM)
set.seed(20260817)
gsL <- intersect(PANEL, names(EX$LIHC)); gsU <- intersect(PANEL, names(EX$LUAD))
gs_common <- intersect(gsL, gsU)
xL <- EX$LIHC$CYB5R3; xU <- EX$LUAD$CYB5R3
ML <- as.matrix(EX$LIHC[, gs_common, drop = FALSE])
MU <- as.matrix(EX$LUAD[, gs_common, drop = FALSE])
rL <- apply(ML, 2, rank); rU <- apply(MU, 2, rank)
for (b in seq_len(N_PERM)) {
  rl <- cor(rank(sample(xL)), rL, use = "pairwise.complete.obs")[1, ]
  ru <- cor(rank(sample(xU)), rU, use = "pairwise.complete.obs")[1, ]
  perm_rev[b] <- sum(sign(rl) != sign(ru), na.rm = TRUE)
}
p_perm <- (sum(perm_rev >= obs_rev) + 1) / (N_PERM + 1)
PERM <- data.frame(n_genes = length(gs_common), observed_reversals = obs_rev,
                   perm_mean = mean(perm_rev), perm_sd = sd(perm_rev),
                   perm_q95 = unname(quantile(perm_rev, .95)),
                   n_permutations = N_PERM, p_permutation = p_perm)
print(PERM); w_res(PERM, "12_permutation_reversals.csv")
message(sprintf("\n  관측 역전 %d개 / 순열 평균 %.1f개 → permutation p = %.4f",
                obs_rev, mean(perm_rev), p_perm))
message("  해석 지침: p_permutation 이 유의하지 않으면 '조직 역전'은 우연과 구별되지 않으므로")
message("  원고에서 '모델과 일치(consistent with)' 수준으로만 서술해야 합니다 (심사 R3 지시).")

## =====================================================================
## 7. 그림
## =====================================================================
try({
  pl <- DIF
  pl$gene <- factor(pl$gene, levels = pl$gene[order(pl$rho_LIHC)])
  long <- rbind(data.frame(gene = pl$gene, cohort = "LIHC", rho = pl$rho_LIHC,
                           lo = pl$lo_LIHC, hi = pl$hi_LIHC),
                data.frame(gene = pl$gene, cohort = "LUAD", rho = pl$rho_LUAD,
                           lo = pl$lo_LUAD, hi = pl$hi_LUAD))
  g <- ggplot2::ggplot(long, ggplot2::aes(rho, gene, colour = cohort)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi),
                            height = .25, position = ggplot2::position_dodge(.5)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(.5), size = 2) +
    ggplot2::scale_colour_manual(values = c(LIHC = "#C0392B", LUAD = "#2E6DA4")) +
    ggplot2::labs(x = "Spearman correlation with CYB5R3 (95% CI)", y = NULL,
                  title = "Exploratory pairwise associations of CYB5R3, liver versus lung") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(FIG, "12_coexpression_liver_lung.png"), g,
                  width = 7.5, height = 6, dpi = 300)
}, silent = TRUE)

save_session("12_coexpression")
message("\n[12] 완료.")
