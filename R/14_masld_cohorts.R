# =====================================================================
#  14_masld_cohorts.R    ★ 본 개정의 핵심 신규 분석 ★
#  목적 (심사 대응):
#   · R3 "Either add fatty-liver data or remove that claim from the title and
#     abstract … Public fibrosis-staged liver datasets exist and would let you
#     test CYB5R3, MTARC1, SCD1 and NQO1 across disease severity … This is the
#     single most valuable addition available to you."
#   · R1-2 "The MASLD–HCC stage-switch model is not directly tested."
#   · R2 "The analyzed resources do not include samples specific to the
#     inflammatory/metabolic form of liver disease."
#
#  코호트 (모두 사람 간생검, 병리의사 판독 조직학 등급 포함)
#   · GSE130970  Hoang SA et al. Sci Rep 2019;9:12541.            n = 78
#   · GSE135251  Govaere O et al. Sci Transl Med 2020;12:eaba4448. n = 216
#   · GSE162694  Pantano L et al. Sci Rep 2021;11:18045.           n = 143
#
#  검정: 조직학 등급(순서형) 대비 Spearman 경향성 + Kruskal–Wallis,
#        3개 코호트 전체에 대해 BH FDR 보정,
#        Fisher-z 랜덤효과 메타분석으로 코호트 간 통합 추정치 산출.
#
#  산출: results/14_masld_trend_all.csv        코호트 x 유전자 x 조직학
#        results/14_masld_meta.csv             메타분석 통합 rho
#        results/14_masld_partial.csv          다변량/편상관 (지방증·풍선변성·염증 분리)
#        results/14_masld_sample_level.csv     환자 수준 병합자료 (원자료 공개용)
#        figures/14_masld_severity.png
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
hdr("14 · 사람 MASLD 간생검 코호트 3종")

options(timeout = 3600)
GEO_DIR <- file.path(CACHE, "geo"); dir.create(GEO_DIR, showWarnings = FALSE, recursive = TRUE)
## 이미 내려받아 둔 자료가 있으면 그대로 씁니다 (재다운로드 방지).
if (dir.exists(EXTERNAL_GEO_CACHE)) {
  message("  기존 GEO 캐시를 사용합니다: ", EXTERNAL_GEO_CACHE)
  GEO_DIR <- EXTERNAL_GEO_CACHE
}

## 유전자 식별자 사전 (심볼 / Entrez / Ensembl 모두 대응)
IDMAP <- list(
  CYB5R3 = c("CYB5R3","DIA1","1727","ENSG00000100243"),
  MTARC1 = c("MTARC1","MARC1","MOSC1","64757","ENSG00000186205"),
  SCD    = c("SCD","SCD1","6319","ENSG00000099194"),
  NQO1   = c("NQO1","1728","ENSG00000181019"),
  PARP16 = c("PARP16","54956","ENSG00000138617"),
  FASN   = c("FASN","2194","ENSG00000169710"),
  COL1A1 = c("COL1A1","1277","ENSG00000108821"),
  TGFB1  = c("TGFB1","7040","ENSG00000105329"))

## ---------------------------------------------------------------------
## 헬퍼
## ---------------------------------------------------------------------
get_pheno <- function(gse) {
  f <- file.path(GEO_DIR, paste0(gse, "_pheno.rds"))
  if (file.exists(f)) return(readRDS(f))
  e <- GEOquery::getGEO(gse, GSEMatrix = TRUE, getGPL = FALSE, destdir = GEO_DIR)
  ph <- Biobase::pData(e[[1]]); saveRDS(ph, f); ph
}
parse_chars <- function(ph) {
  cc <- grep("^characteristics_ch1", names(ph), value = TRUE)
  out <- ph["geo_accession"]; out$title <- as.character(ph$title)
  for (cn in cc) {
    v <- as.character(ph[[cn]])
    key <- unique(trimws(sub(":.*$", "", v))); key <- key[nzchar(key)][1]
    if (is.na(key)) next
    out[[make.names(key)]] <- trimws(sub("^[^:]*:\\s*", "", v))
  }
  out
}
supp_url <- function(gse) sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/suppl/",
                                  substr(gse, 1, nchar(gse) - 3), gse)

#' 발현행렬 확보: (1) 캐시된 단일 매트릭스 → (2) FTP 보충파일 → (3) per-sample counts 폴더
load_matrix <- function(gse) {
  d <- file.path(GEO_DIR, gse); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  pat <- "\\.(csv|tsv|txt)(\\.gz)?$"
  have <- c(list.files(d, full.names = TRUE, pattern = pat),
            list.files(GEO_DIR, full.names = TRUE, pattern = paste0("^", gse, ".*", pat)))
  have <- have[!grepl("series_matrix", have)]
  cand <- have[grepl("tpm|fpkm|count|expr|matrix|salmon", have, ignore.case = TRUE)]
  if (!length(cand)) {
    url <- supp_url(gse)
    html <- tryCatch(paste(readLines(url, warn = FALSE), collapse = "\n"),
                     error = function(e) NA_character_)
    fn <- if (!is.na(html)) unique(regmatches(html,
            gregexpr("GSE\\d+_[^\"'> ]+\\.(csv|tsv|txt)\\.gz", html, perl = TRUE))[[1]]) else character(0)
    for (f in fn) {
      dest <- file.path(d, f)
      if (!file.exists(dest)) tryCatch(download.file(paste0(url, f), dest, mode = "wb", quiet = TRUE),
                                       error = function(e) NULL)
    }
    cand <- list.files(d, full.names = TRUE, pattern = pat)
    cand <- cand[grepl("tpm|fpkm|count|expr|matrix|salmon", cand, ignore.case = TRUE)]
  }
  if (length(cand)) {
    norm <- cand[grepl("tpm|fpkm", cand, ignore.case = TRUE)]
    pick <- if (length(norm)) norm[1] else cand[1]
    m <- as.data.frame(data.table::fread(pick, check.names = FALSE))
    attr(m, "is_norm") <- length(norm) > 0
    attr(m, "source")  <- basename(pick)
    message("  행렬: ", basename(pick), "  (", nrow(m), " x ", ncol(m), ")")
    return(m)
  }
  ## per-sample counts 폴더 (예: GSE135251/counts/GSM*.counts.txt.gz)
  cdir <- list.dirs(d, recursive = TRUE)
  cdir <- cdir[grepl("count", cdir, ignore.case = TRUE)]
  fs <- if (length(cdir)) list.files(cdir[1], full.names = TRUE, pattern = "\\.txt(\\.gz)?$") else character(0)
  if (length(fs) > 20) {
    message("  per-sample counts ", length(fs), "개를 병합합니다 …")
    lst <- lapply(fs, function(f) {
      x <- data.table::fread(f, header = FALSE)
      data.table::setnames(x, c("id","count"))
      x$sample <- sub("_.*$", "", basename(f)); x
    })
    all <- data.table::rbindlist(lst)
    m <- data.table::dcast(all, id ~ sample, value.var = "count")
    m <- as.data.frame(m); attr(m, "is_norm") <- FALSE
    attr(m, "source") <- paste0(basename(cdir[1]), " (", length(fs), " files merged)")
    return(m)
  }
  message("  [실패] ", gse, " 발현행렬을 찾지 못했습니다. 아래에서 받아 ", d, " 에 넣으십시오:\n    ", supp_url(gse))
  NULL
}
to_cpm <- function(m) {
  if (is.null(m)) return(NULL)
  if (isTRUE(attr(m, "is_norm"))) return(m)
  v <- as.matrix(m[, -1, drop = FALSE]); mode(v) <- "numeric"
  v <- sweep(v, 2, colSums(v, na.rm = TRUE), "/") * 1e6
  out <- cbind(m[, 1, drop = FALSE], as.data.frame(v))
  attr(out, "is_norm") <- TRUE; attr(out, "source") <- attr(m, "source"); out
}
extract_genes <- function(m, genes = names(IDMAP)) {
  if (is.null(m)) return(NULL)
  id <- toupper(sub("\\..*$", "", as.character(m[[1]])))
  out <- data.frame(sample = colnames(m)[-1])
  for (g in genes) {
    i <- which(id %in% toupper(IDMAP[[g]]))
    if (!length(i)) { message("    [miss] ", g); next }
    v <- suppressWarnings(as.numeric(m[i[1], -1]))
    if (all(is.na(v))) next
    out[[g]] <- log2(v + 1)
  }
  out
}
trend_test <- function(x, grade) {
  ok <- !is.na(x) & !is.na(grade); x <- x[ok]; g <- grade[ok]
  if (length(unique(g)) < 2 || length(x) < 10) return(NULL)
  sp <- suppressWarnings(cor.test(x, as.numeric(g), method = "spearman", exact = FALSE))
  kw <- kruskal.test(x, factor(g))
  data.frame(n = length(x), n_groups = length(unique(g)),
             rho = unname(sp$estimate), p_trend = sp$p.value, p_kruskal = kw$p.value)
}
## 조직학 변수명 정규화 — 코호트마다 표기가 달라 공통 축으로 사상한다.
canon_histology <- function(nm) {
  n <- tolower(nm)
  dplyr::case_when(
    grepl("fibrosis", n)                              ~ "fibrosis stage",
    grepl("nas|activity", n)                          ~ "NAFLD activity score",
    grepl("steato", n)                                ~ "steatosis grade",
    grepl("balloon", n)                               ~ "ballooning grade",
    grepl("lobular|inflamm", n)                       ~ "lobular inflammation grade",
    TRUE ~ NA_character_)
}

## ---------------------------------------------------------------------
## 코호트별 처리
## ---------------------------------------------------------------------
COHORTS <- c("GSE130970","GSE135251","GSE162694")
trend_rows <- list(); sample_rows <- list(); partial_rows <- list()

for (gse in COHORTS) {
  hdr(paste("코호트", gse))
  m <- tryCatch(to_cpm(load_matrix(gse)), error = function(e) { message("  ", conditionMessage(e)); NULL })
  if (is.null(m)) next
  ex <- extract_genes(m)
  if (is.null(ex) || !"CYB5R3" %in% names(ex)) { message("  CYB5R3 추출 실패 — 건너뜀"); next }

  ph <- get_pheno(gse); pc <- parse_chars(ph)
  ## ---- 샘플명 매칭: 후보 표의 모든 열을 정규화해 가장 잘 맞는 열을 고른다 ----
  norm_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))
  cand_tabs <- list(characteristics = pc, pData = as.data.frame(ph))
  alt <- c(list.files(GEO_DIR, pattern = paste0("^", gse, ".*(pheno|meta|clin).*\\.csv$"),
                      full.names = TRUE),
           list.files(file.path(GEO_DIR, gse), pattern = "(pheno|meta|clin).*\\.csv$",
                      full.names = TRUE))
  for (a in alt) cand_tabs[[basename(a)]] <- as.data.frame(data.table::fread(a))
  exn <- norm_id(ex$sample)
  best <- list(frac = 0)
  for (tn in names(cand_tabs)) {
    tb <- cand_tabs[[tn]]
    for (k in names(tb)) {
      fr <- mean(exn %in% norm_id(tb[[k]]))
      if (fr > best$frac) best <- list(frac = fr, table = tn, col = k)
    }
  }
  ## 2차: 정확 일치가 낮으면 부분문자열(포함관계)로 재시도한다.
  ##      GSE162694 처럼 발현행렬 열이 '548nash1' 이고 pheno 의 title 이
  ##      그 문자열을 포함하거나 그 안에 포함되는 경우를 잡는다.
  if (best$frac < 0.5) {
    for (tn in names(cand_tabs)) {
      tb <- cand_tabs[[tn]]
      for (k in names(tb)) {
        vals <- norm_id(tb[[k]])
        if (!any(nzchar(vals))) next
        hit <- vapply(exn, function(e) {
          if (!nzchar(e)) return(FALSE)
          w <- which(vals == e | grepl(e, vals, fixed = TRUE) | 
                     (nchar(vals) > 2 & vapply(vals, function(v) grepl(v, e, fixed = TRUE), TRUE)))
          length(w) == 1L
        }, TRUE)
        fr <- mean(hit)
        if (fr > best$frac) best <- list(frac = fr, table = tn, col = k, fuzzy = TRUE)
      }
    }
  }
  message(sprintf("  [match] 최적 키: %s$%s  (일치율 %.1f%%%s)",
                  best$table %||% "-", best$col %||% "-", 100 * best$frac,
                  if (isTRUE(best$fuzzy)) ", 부분일치" else ""))
  d <- data.frame()
  if (best$frac >= 0.5 && isTRUE(best$fuzzy)) {
    tb <- cand_tabs[[best$table]]
    vals <- norm_id(tb[[best$col]])
    idx <- vapply(exn, function(e) {
      w <- which(vals == e | grepl(e, vals, fixed = TRUE) |
                 (nchar(vals) > 2 & vapply(vals, function(v) grepl(v, e, fixed = TRUE), TRUE)))
      if (length(w) == 1L) w else NA_integer_
    }, NA_integer_)
    keep <- !is.na(idx)
    d <- cbind(ex[keep, , drop = FALSE], tb[idx[keep], , drop = FALSE])
    if (best$table != "characteristics" && "geo_accession" %in% names(d))
      d <- merge(d, pc, by = "geo_accession", all.x = TRUE, suffixes = c("", ".pc"))
  } else if (best$frac >= 0.5) {
    tb <- cand_tabs[[best$table]]
    tb$.key <- norm_id(tb[[best$col]])
    ex2 <- ex; ex2$.key <- exn
    d <- merge(ex2, tb, by = ".key")
    ## 특성표(pc)가 아닌 표로 매칭됐다면 조직학 열이 없을 수 있으므로 pc 를 덧붙인다
    if (best$table != "characteristics" && "geo_accession" %in% names(d))
      d <- merge(d, pc, by = "geo_accession", all.x = TRUE, suffixes = c("", ".pc"))
  }
  message("  매칭 표본: ", nrow(d))
  if (!nrow(d)) {
    message("  [진단] 발현행렬 열 이름 예시: ", paste(utils::head(ex$sample, 5), collapse = ", "))
    for (tn in names(cand_tabs))
      message("  [진단] ", tn, " 열: ", paste(utils::head(names(cand_tabs[[tn]]), 12), collapse = ", "))
    next
  }

  hcols <- grep("fibrosis|nas|activity|steato|balloon|lobular|inflamm", names(d),
                ignore.case = TRUE, value = TRUE)
  hcols <- hcols[!grepl("^CYB5R3$|^MTARC1$", hcols)]
  message("  조직학 변수: ", paste(hcols, collapse = ", "))
  genes_here <- intersect(names(IDMAP), names(d))

  for (h in hcols) {
    gv <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(d[[h]]))))
    if (all(is.na(gv)) || length(unique(na.omit(gv))) < 2) next
    can <- canon_histology(h); if (is.na(can)) next
    for (g in genes_here) {
      t <- trend_test(d[[g]], gv); if (is.null(t)) next
      trend_rows[[length(trend_rows)+1]] <-
        cbind(data.frame(cohort = gse, gene = g, histology_raw = h, histology = can), t)
    }
    d[[paste0("h_", make.names(can))]] <- gv
  }

  ## 조직학 성분 분리: 지방증 · 소엽염증 · 풍선변성을 동시에 넣은 선형모형
  need <- c("h_steatosis.grade","h_lobular.inflammation.grade","h_ballooning.grade")
  if (all(need %in% names(d))) {
    for (g in genes_here) {
      fit <- lm(as.formula(sprintf("%s ~ h_steatosis.grade + h_lobular.inflammation.grade + h_ballooning.grade", g)),
                data = d)
      s <- summary(fit)$coefficients
      for (tm in rownames(s)[-1])
        partial_rows[[length(partial_rows)+1]] <- data.frame(
          cohort = gse, gene = g, model = "expr ~ steatosis + lobular inflammation + ballooning",
          term = tm, beta = s[tm,1], se = s[tm,2], p = s[tm,4])
      ## 편상관 (다른 두 성분 보정)
      for (tm in need) {
        others <- setdiff(need, tm)
        ps <- partial_spearman(d[[g]], d[[tm]], d[, others, drop = FALSE])
        partial_rows[[length(partial_rows)+1]] <- data.frame(
          cohort = gse, gene = g, model = paste("partial Spearman adjusting for", paste(others, collapse = " + ")),
          term = tm, beta = ps$rho, se = NA_real_, p = ps$p)
      }
    }
  }

  keep <- c("sample", genes_here, grep("^h_", names(d), value = TRUE))
  sample_rows[[gse]] <- cbind(cohort = gse, d[, keep, drop = FALSE])
}

## ---------------------------------------------------------------------
## 통합 · FDR · 메타분석
## ---------------------------------------------------------------------
TR <- dplyr::bind_rows(trend_rows)
if (!nrow(TR)) stop("경향성 결과가 없습니다 — 위 로그에서 자료 인식 상태를 확인하십시오.")
TR$q_trend <- p.adjust(TR$p_trend, method = "BH")   # 코호트·유전자·조직학 전체에 대해 보정
TR <- TR[order(TR$gene, TR$histology, TR$cohort), ]
print(as.data.frame(TR), digits = 3); w_res(TR, "14_masld_trend_all.csv")

## Fisher-z 랜덤효과 메타분석 (유전자 x 조직학)
meta_rows <- list()
for (g in unique(TR$gene)) for (h in unique(TR$histology)) {
  s <- TR[TR$gene == g & TR$histology == h, ]
  s <- s[!is.na(s$rho) & s$n > 5, ]
  if (nrow(s) < 2) next
  z <- atanh(pmin(pmax(s$rho, -0.999), 0.999)); vi <- 1 / (s$n - 3)
  est <- if (has_pkg("metafor")) {
    fit <- metafor::rma(yi = z, vi = vi, method = "REML")
    data.frame(rho_pooled = tanh(as.numeric(fit$b)),
               lo = tanh(fit$ci.lb), hi = tanh(fit$ci.ub),
               p = fit$pval, I2 = fit$I2, tau2 = fit$tau2)
  } else {
    w <- 1 / vi; mu <- sum(w * z) / sum(w); se <- sqrt(1 / sum(w))
    data.frame(rho_pooled = tanh(mu), lo = tanh(mu - 1.96 * se), hi = tanh(mu + 1.96 * se),
               p = 2 * pnorm(-abs(mu / se)), I2 = NA, tau2 = NA)
  }
  meta_rows[[length(meta_rows)+1]] <- cbind(
    data.frame(gene = g, histology = h, k_cohorts = nrow(s), n_total = sum(s$n)), est)
}
MT <- dplyr::bind_rows(meta_rows)
if (nrow(MT)) { MT$q <- p.adjust(MT$p, method = "BH"); print(as.data.frame(MT), digits = 3)
                w_res(MT, "14_masld_meta.csv") }

if (length(partial_rows)) {
  PR <- dplyr::bind_rows(partial_rows); PR$q <- p.adjust(PR$p, method = "BH")
  print(as.data.frame(PR), digits = 3); w_res(PR, "14_masld_partial.csv")
}
SR <- dplyr::bind_rows(sample_rows)
if (nrow(SR)) w_res(SR, "14_masld_sample_level.csv")

## ---------------------------------------------------------------------
## 그림 — 조직학 중증도에 따른 발현
## ---------------------------------------------------------------------
try({
  pl <- MT[MT$gene %in% c("CYB5R3","MTARC1","SCD","NQO1"), ]
  pl$gene <- factor(pl$gene, levels = c("CYB5R3","MTARC1","SCD","NQO1"))
  g <- ggplot2::ggplot(pl, ggplot2::aes(rho_pooled, histology, colour = gene)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), height = .2,
                            position = ggplot2::position_dodge(.6)) +
    ggplot2::geom_point(size = 2.4, position = ggplot2::position_dodge(.6)) +
    ggplot2::labs(x = "Pooled Spearman correlation with histological grade (95% CI)",
                  y = NULL, colour = NULL,
                  title = "CYB5R3, MTARC1, SCD and NQO1 across the human MASLD spectrum",
                  subtitle = sprintf("Random-effects meta-analysis of %d liver-biopsy cohorts",
                                     dplyr::n_distinct(TR$cohort))) +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(FIG, "14_masld_severity.png"), g, width = 8, height = 5, dpi = 300)
}, silent = TRUE)

save_session("14_masld")
message("\n[14] 완료.")
