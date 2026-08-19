# =====================================================================
#  15_depmap.R
#  목적 (심사 대응):
#   · R3 "The DepMap comparison is p = 0.17 in the main text and p = 0.165 in
#     the supplement." → 하나의 산출값을 만들고 본문·보충 모두 여기서 인용한다.
#   · R1-3 "the DepMap analysis does not identify CYB5R3 as a strong genetic
#     dependency" → 음성 결과를 Results 본문에 명시적으로 보고하기 위한 수치.
#  산출: results/15_depmap_dependency.csv
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
hdr("15 · DepMap CRISPR 의존성 (liver vs lung)")

get_depmap <- function() {
  ## (A) Bioconductor depmap (버전 고정)
  if (requireNamespace("depmap", quietly = TRUE) &&
      requireNamespace("ExperimentHub", quietly = TRUE)) {
    ok <- tryCatch({
      eh <- ExperimentHub::ExperimentHub()
      ce <- depmap::depmap_crispr()
      md <- depmap::depmap_metadata()
      list(crispr = ce, meta = md,
           source = sprintf("Bioconductor depmap %s",
                            as.character(utils::packageVersion("depmap"))))
    }, error = function(e) { message("  depmap 패키지 실패: ", conditionMessage(e)); NULL })
    if (!is.null(ok)) return(ok)
  }
  ## (B) 수동 경로: DepMap 포털에서 받은 CSV
  f1 <- file.path(CACHE, "CRISPRGeneEffect.csv"); f2 <- file.path(CACHE, "Model.csv")
  if (file.exists(f1) && file.exists(f2)) {
    ge <- as.data.frame(data.table::fread(f1))
    md <- as.data.frame(data.table::fread(f2))
    col <- grep("^CYB5R3\\b|^CYB5R3 \\(", names(ge), value = TRUE)[1]
    ce <- data.frame(depmap_id = ge[[1]], dependency = suppressWarnings(as.numeric(ge[[col]])))
    lin <- grep("OncotreeLineage|^lineage$", names(md), value = TRUE)[1]
    md2 <- data.frame(depmap_id = md[[grep("ModelID|depmap_id", names(md), value = TRUE)[1]]],
                      lineage = md[[lin]])
    return(list(crispr = ce, meta = md2, source = "DepMap portal CSV (CRISPRGeneEffect.csv + Model.csv)"))
  }
  stop("DepMap 자료를 찾지 못했습니다.\n",
       "  BiocManager::install('depmap') 하시거나,\n",
       "  DepMap 포털에서 CRISPRGeneEffect.csv 와 Model.csv 를 받아 ", CACHE, " 에 넣으십시오.")
}

dm <- get_depmap()
ce <- dm$crispr; md <- dm$meta
if ("gene_name" %in% names(ce)) ce <- ce[ce$gene_name == "CYB5R3", ]
key <- intersect(names(ce), names(md)); key <- key[key %in% c("depmap_id")][1]
D <- merge(ce, md, by = key)
lincol <- grep("lineage", names(D), value = TRUE)[1]
depcol <- grep("^dependency$|gene_effect", names(D), value = TRUE)[1]
D$grp <- dplyr::case_when(
  grepl("liver", D[[lincol]], ignore.case = TRUE) ~ "Liver (HCC)",
  grepl("lung",  D[[lincol]], ignore.case = TRUE) ~ "Lung",
  TRUE ~ NA_character_)
D <- D[!is.na(D$grp) & !is.na(D[[depcol]]), ]
message("  lines: ", paste(names(table(D$grp)), table(D$grp), collapse = " | "))

w <- wilcox.test(D[[depcol]] ~ D$grp)
out <- D %>% dplyr::group_by(grp) %>% dplyr::summarise(
  n_lines = dplyr::n(),
  median_gene_effect = median(.data[[depcol]]),
  Q1 = quantile(.data[[depcol]], .25), Q3 = quantile(.data[[depcol]], .75),
  mean_gene_effect = mean(.data[[depcol]]),
  pct_below_minus0.5 = 100 * mean(.data[[depcol]] < -0.5), .groups = "drop")
out$wilcoxon_p <- signif(w$p.value, 3)
out$median_difference_liver_minus_lung <-
  out$median_gene_effect[out$grp == "Liver (HCC)"] - out$median_gene_effect[out$grp == "Lung"]
out$data_source <- dm$source
out$note <- "Chronos gene effect: 0 = no effect, -1 = median common-essential, < -0.5 = dependent line"
print(as.data.frame(out), digits = 4); w_res(out, "15_depmap_dependency.csv")
message(sprintf("\n  ★ 본문·보충 모두 이 p 값을 인용하십시오: Wilcoxon p = %s", signif(w$p.value, 3)))

save_session("15_depmap")
message("\n[15] 완료.")
