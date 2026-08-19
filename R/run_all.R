# =====================================================================
#  run_all.R — 전체 파이프라인 실행
#  사용법:  RStudio 에서 이 파일을 열고 Source,  또는 터미널에서
#           Rscript R/run_all.R
#  소요시간 참고: 12번(순열 10,000회)과 14번(GEO 다운로드)이 가장 오래 걸립니다.
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
R_HOME_DIR <- .here()
source(file.path(R_HOME_DIR, "01_common.R"))

STEPS <- c("10_expression_rawdata.R",
           "11_survival_primary.R",
           "12_coexpression_tissue.R",
           "13_pancancer.R",
           "14_masld_cohorts.R",
           "15_depmap.R",
           "16_lihc_stratified.R",
           "17_manuscript_figures.R")

## 특정 단계만 돌리려면 예: RUN_ONLY <- c("14_masld_cohorts.R")
if (!exists("RUN_ONLY")) RUN_ONLY <- STEPS

log_file <- file.path(RES, sprintf("run_all_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
con <- file(log_file, open = "wt"); sink(con, split = TRUE); sink(con, type = "message")
on.exit({ sink(type = "message"); sink(); close(con) }, add = TRUE)

status <- list()
for (s in intersect(STEPS, RUN_ONLY)) {
  t0 <- Sys.time()
  message("\n\n############ ", s, " ############")
  ok <- tryCatch({ source(file.path(R_HOME_DIR, s), local = new.env()); TRUE },
                 error = function(e) { message("!! 실패: ", conditionMessage(e)); FALSE })
  status[[s]] <- data.frame(script = s, ok = ok,
                            minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2))
}
st <- dplyr::bind_rows(status); print(st)
data.table::fwrite(st, file.path(RES, "run_all_status.csv"))
writeLines(capture.output(sessionInfo()), file.path(RES, "sessionInfo_run_all.txt"))
message("\n로그: ", log_file)
