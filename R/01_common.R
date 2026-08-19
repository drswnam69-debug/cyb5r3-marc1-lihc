# =====================================================================
#  CYB5R3–mARC1 / IJMS ijms-4481300 — Revision 1 reproducible pipeline
#  01_common.R : 경로 · 패키지 · 공통 헬퍼 (직접 실행하지 말고 source() 하십시오)
#
#  이 파일 하나만 고치면 모든 스크립트의 경로/설정이 바뀝니다.
#  작성: 2026-08  (심사 R1 대응)
# =====================================================================

## ---- 0. 재현성 고정 ------------------------------------------------
set.seed(20260817)
options(stringsAsFactors = FALSE, scipen = 999)

## ---- 1. 패키지 -----------------------------------------------------
.cran <- c("data.table","dplyr","tidyr","stringr","readr","tibble","purrr",
           "ggplot2","ggpubr","survival","survminer","boot","R.utils")
.cran_opt <- c("rms","mice","metafor","timeROC","pheatmap","forcats","ggrepel")
.bioc <- c("UCSCXenaTools","GEOquery","Biobase")

.load <- function(pkgs, quiet = TRUE) {
  ok <- vapply(pkgs, function(p) {
    if (requireNamespace(p, quietly = TRUE)) { suppressPackageStartupMessages(
      library(p, character.only = TRUE)); TRUE } else FALSE
  }, logical(1))
  if (!quiet) for (i in seq_along(pkgs))
    message(sprintf("  %-16s %s", pkgs[i], if (ok[i]) "OK" else "MISSING"))
  invisible(ok)
}
.load(c(.cran, .bioc))
.OPT <- .load(.cran_opt)          # 선택 패키지: 없으면 해당 분석만 건너뜀
has_pkg <- function(p) isTRUE(unname(.OPT[p])) || requireNamespace(p, quietly = TRUE)

## ---- 2. 경로 (스크립트 위치 자동 탐지) -----------------------------
.this_dir <- function() {
  fr <- sys.frames()
  if (length(fr)) for (i in rev(seq_along(fr))) {
    o <- fr[[i]]$ofile
    if (!is.null(o) && nzchar(o)) return(dirname(normalizePath(o)))
  }
  m <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(m)) return(dirname(normalizePath(sub("^--file=", "", m[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- try(rstudioapi::getSourceEditorContext()$path, silent = TRUE)
    if (!inherits(p, "try-error") && !is.null(p) && nzchar(p))
      return(dirname(normalizePath(p)))
  }
  getwd()
}
R_DIR <- .this_dir()
ROOT  <- normalizePath(file.path(R_DIR, ".."), mustWork = FALSE)
RES   <- file.path(ROOT, "results")
FIG   <- file.path(ROOT, "figures")
CACHE <- file.path(ROOT, "cache")
for (d in c(RES, FIG, CACHE)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

## 이미 내려받아 둔 GEO 원자료가 있으면 재사용 (없으면 새로 내려받음).
## 필요하면 아래 경로를 교수님 환경에 맞게 바꾸십시오.
EXTERNAL_GEO_CACHE <- Sys.getenv(
  "CYB5R3_GEO_CACHE",
  unset = normalizePath(file.path(ROOT, "..", "Revision_2026_human_MASLD",
                                  "analysis", "cache", "geo"), mustWork = FALSE))

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

## ---- 3. 유전자 정의 -------------------------------------------------
## 명명법 (심사 minor 7): 유전자 = CYB5R3 / MTARC1 / SCD, 단백질 = CYB5R3 / mARC1 / SCD1
AXIS   <- c("CYB5R3", "MTARC1")
ALIAS  <- list(CYB5R3 = c("CYB5R3","DIA1"),
               MTARC1 = c("MTARC1","MARC1","MOSC1"),
               MTARC2 = c("MTARC2","MARC2","MOSC2"),
               SCD    = c("SCD","SCD1"))

## 사전 지정 16-유전자 redox/antioxidant 패널 (원 투고 Supplementary S5 와 동일)
REDOX_PANEL <- c("SCD","PARP16","NQO1","NFE2L2","KEAP1","GCLC","GCLM","GSR",
                 "GPX4","TXN","TXNRD1","SOD1","SOD2","CAT","PRDX1","G6PD")

## 종양순도 대리지표용 기질·면역 표지자 (ABSOLUTE 순도를 못 얻을 때만 사용)
STROMAL_IMMUNE <- c("COL1A1","COL1A2","COL3A1","FN1","ACTA2","PDGFRB","THY1","LUM",
                    "PTPRC","CD3E","CD2","CD8A","CD68","CSF1R","ITGAM","LYZ","HLA-DRA")

## MASLD 코호트에서 볼 유전자 (심사 R3: CYB5R3, MTARC1, SCD1, NQO1 across severity)
MASLD_GENES <- c("CYB5R3","MTARC1","SCD","NQO1","PARP16","FASN","COL1A1","TGFB1")

## ---- 4. Xena 호스트/데이터셋 ----------------------------------------
XH <- list(
  tcga = "https://tcga.xenahubs.net",
  toil = "https://toil.xenahubs.net",
  pan  = "https://pancanatlas.xenahubs.net"
)
XDS <- list(
  lihc_hiseq  = "TCGA.LIHC.sampleMap/HiSeqV2",
  luad_hiseq  = "TCGA.LUAD.sampleMap/HiSeqV2",
  lihc_clin   = "TCGA.LIHC.sampleMap/LIHC_clinicalMatrix",
  toil_tpm    = "TcgaTargetGtex_rsem_gene_tpm",
  toil_pheno  = "TcgaTargetGTEX_phenotype.txt",
  pan_expr    = "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena",
  pan_surv    = "Survival_SupplementalTable_S1_20171025_xena_sp",
  pan_purity  = "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"
)

## ---- 5. 바코드 헬퍼 -------------------------------------------------
nb   <- function(x) toupper(gsub("\\.", "-", as.character(x)))
pat  <- function(x) substr(nb(x), 1, 12)          # 환자 12자리
typ  <- function(x) substr(nb(x), 14, 15)         # 샘플 유형 2자리
bc15 <- function(x) substr(nb(x), 1, 15)
TUMOR_CODES  <- c("01","02","03","05","06","07")
NORMAL_CODES <- c("11")

tovec <- function(m) {
  if (is.null(m)) return(NULL)
  if (is.null(dim(m))) return(m)
  v <- as.numeric(m[1, ]); names(v) <- colnames(m); v
}

#' Xena 단일 유전자 조회 (별칭 순차 시도). 실패 시 NULL.
fetch_gene <- function(host, ds, gene, probemap = TRUE) {
  tries <- ALIAS[[gene]] %||% gene
  for (t in tries) {
    m <- tryCatch(UCSCXenaTools::fetch_dense_values(host, ds, t,
                                                    use_probeMap = probemap),
                  error = function(e) NULL)
    if (!is.null(m) && !all(is.na(m))) {
      attr(m, "alias_used") <- t
      return(m)
    }
  }
  NULL
}

#' Xena 데이터셋을 '표'로 내려받기 (범주형 열 보존).
#' 정확한 ID 를 모를 때를 대비해 XenaData 에서 정규식으로 후보를 찾아 진단 출력한다.
xena_table <- function(pattern, exact = NULL, hub = NULL) {
  xd <- UCSCXenaTools::XenaData
  if (!is.null(hub) && "XenaHostNames" %in% names(xd))
    xd <- xd[xd$XenaHostNames %in% hub, , drop = FALSE]
  cand <- unique(xd$XenaDatasets[grepl(pattern, xd$XenaDatasets, ignore.case = TRUE)])
  if (!is.null(exact) && exact %in% xd$XenaDatasets) cand <- c(exact, setdiff(cand, exact))
  if (!length(cand)) {
    message("  [xena_table] '", pattern, "' 에 맞는 데이터셋이 없습니다.")
    return(NULL)
  }
  message("  [xena_table] 후보 ", length(cand), "개; 사용: ", cand[1])
  if (length(cand) > 1) message("     (기타: ", paste(utils::head(cand[-1], 4), collapse = " | "), ")")
  out <- tryCatch({
    q <- UCSCXenaTools::XenaGenerate(subset = XenaDatasets == cand[1])
    f <- UCSCXenaTools::XenaQuery(q) |>
         UCSCXenaTools::XenaDownload(destdir = CACHE, trans_slash = TRUE, force = FALSE)
    p <- UCSCXenaTools::XenaPrepare(f)
    if (is.list(p) && !is.data.frame(p)) p <- p[[1]]
    as.data.frame(p)
  }, error = function(e) { message("  [xena_table] 실패: ", conditionMessage(e)); NULL })
  out
}

#' 2수준 이상일 때만 factor 를 돌려주고, 아니면 NULL (모형에서 제외)
safe_factor <- function(x) {
  x <- as.character(x); x[!nzchar(trimws(x))] <- NA
  if (length(unique(na.omit(x))) < 2) return(NULL)
  factor(x)
}

#' 여러 유전자 → sample x gene data.frame
fetch_matrix <- function(host, ds, genes, probemap = TRUE) {
  out <- list(); used <- character(0)
  for (g in genes) {
    m <- fetch_gene(host, ds, g, probemap)
    if (is.null(m)) { message("  [skip] ", g, " 조회 실패: ", ds); next }
    out[[g]] <- tovec(m); used[g] <- attr(m, "alias_used")
  }
  if (!length(out)) return(NULL)
  ln <- unique(lengths(out))
  if (length(ln) > 1) {           # 길이가 다르면 이름 기준 정렬
    allnm <- Reduce(union, lapply(out, names))
    out <- lapply(out, function(v) v[allnm])
  }
  df <- as.data.frame(out, check.names = FALSE)
  df$sample <- nb(names(out[[1]]))
  attr(df, "alias_used") <- used
  df
}

## ---- 6. 통계 헬퍼 ---------------------------------------------------
#' Spearman rho + Fisher z 기반 95% CI + p
spearman_ci <- function(x, y, conf = 0.95) {
  ok <- complete.cases(x, y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 5) return(data.frame(n = n, rho = NA, lo = NA, hi = NA, p = NA))
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  r  <- unname(ct$estimate)
  z  <- atanh(pmin(pmax(r, -0.999999), 0.999999))
  se <- 1 / sqrt(n - 3)                       # Fieller 보정 없이 보수적으로 사용
  cv <- qnorm(1 - (1 - conf) / 2)
  data.frame(n = n, rho = r,
             lo = tanh(z - cv * se), hi = tanh(z + cv * se),
             p  = ct$p.value)
}

#' 두 독립 코호트의 상관계수 차이 검정 (Fisher r-to-z)  — 심사 R1-3 / R3-2
fisher_z_diff <- function(r1, n1, r2, n2) {
  z1 <- atanh(pmin(pmax(r1, -0.999999), 0.999999))
  z2 <- atanh(pmin(pmax(r2, -0.999999), 0.999999))
  se <- sqrt(1/(n1 - 3) + 1/(n2 - 3))
  z  <- (z1 - z2) / se
  data.frame(z_diff = z, p_diff = 2 * pnorm(-abs(z)),
             d_lo = tanh((z1 - z2) - qnorm(0.975) * se),
             d_hi = tanh((z1 - z2) + qnorm(0.975) * se))
}

#' 편(partial) Spearman: 순위변환 후 공변량 잔차 상관 — 종양순도 보정용
partial_spearman <- function(x, y, z) {
  z <- as.data.frame(z)
  ok <- complete.cases(x, y, z); x <- x[ok]; y <- y[ok]; z <- z[ok, , drop = FALSE]
  n <- length(x); k <- ncol(z)
  if (n < 10) return(data.frame(n = n, rho = NA, p = NA))
  rx <- residuals(lm(rank(x) ~ ., data = as.data.frame(lapply(z, rank))))
  ry <- residuals(lm(rank(y) ~ ., data = as.data.frame(lapply(z, rank))))
  r  <- cor(rx, ry)
  df <- n - 2 - k
  tt <- r * sqrt(df / (1 - r^2))
  data.frame(n = n, rho = r, p = 2 * pt(-abs(tt), df))
}

#' 역 Kaplan–Meier 로 중앙 추적기간
median_followup <- function(time, event) {
  fit <- survival::survfit(survival::Surv(time, 1 - event) ~ 1)
  s <- summary(fit)$table
  unname(s["median"])
}

#' Harrell C-index + 부트스트랩 ΔC (기저모형 대비)
c_index_delta <- function(data, time, event, base_rhs, full_rhs, B = 1000) {
  f0 <- as.formula(sprintf("Surv(%s,%s) ~ %s", time, event, base_rhs))
  f1 <- as.formula(sprintf("Surv(%s,%s) ~ %s", time, event, full_rhs))
  cc <- function(d) {
    c0 <- summary(survival::coxph(f0, data = d))$concordance[1]
    c1 <- summary(survival::coxph(f1, data = d))$concordance[1]
    c(c0, c1, c1 - c0)
  }
  obs <- cc(data)
  bs <- matrix(NA_real_, B, 3)
  for (b in seq_len(B)) {
    d <- data[sample.int(nrow(data), replace = TRUE), , drop = FALSE]
    bs[b, ] <- tryCatch(cc(d), error = function(e) c(NA, NA, NA))
  }
  q <- function(j) quantile(bs[, j], c(.025, .975), na.rm = TRUE)
  data.frame(C_base = obs[1], C_base_lo = q(1)[1], C_base_hi = q(1)[2],
             C_full = obs[2], C_full_lo = q(2)[1], C_full_hi = q(2)[2],
             dC = obs[3], dC_lo = q(3)[1], dC_hi = q(3)[2], B = B)
}

#' Cliff's delta (효과크기, 비모수)
cliffs_delta <- function(a, b) {
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  if (!length(a) || !length(b)) return(NA_real_)
  ww <- suppressWarnings(wilcox.test(a, b))
  (2 * unname(ww$statistic) / (length(a) * length(b))) - 1
}

#' 요약통계 한 줄 (중앙값·IQR·평균·SD·n)   [주의] dplyr::desc 와 충돌을 피해 이름을 구분
desc_stats <- function(x) {
  x <- x[!is.na(x)]
  data.frame(n = length(x), median = median(x),
             q1 = unname(quantile(x, .25)), q3 = unname(quantile(x, .75)),
             mean = mean(x), sd = sd(x))
}

## ---- 7. 출력 헬퍼 ---------------------------------------------------
w_res <- function(df, name) {
  p <- file.path(RES, name)
  data.table::fwrite(df, p)
  message("  → ", p)
  invisible(p)
}
save_session <- function(tag) {
  writeLines(capture.output(sessionInfo()),
             file.path(RES, sprintf("sessionInfo_%s.txt", tag)))
  writeLines(c(paste("script :", tag),
               paste("run at :", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
               paste("R      :", R.version.string),
               paste("seed   : 20260817")),
             file.path(RES, sprintf("RUNLOG_%s.txt", tag)))
}
hdr <- function(x) message("\n", strrep("=", 68), "\n  ", x, "\n", strrep("=", 68))

message("01_common.R loaded — results → ", RES)
