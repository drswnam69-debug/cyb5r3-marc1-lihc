# =====================================================================
#  16_lihc_stratified.R
#  목적 (심사 대응):
#   · R3 "New results currently appear in the Discussion: move the pan-cancer,
#     grade, aetiology and DepMap analyses into Results, with their methods in
#     Materials and Methods." → 병기·조직등급·병인 분석을 파이프라인에 편입해
#     본문 수치가 전부 재생성되도록 한다.
#   · R3 "Relabel the 'metabolic subset' — negative viral serology is not the
#     same as metabolic disease, and aetiology could be assigned for only about
#     165 of roughly 365 tumours." → serology-negative 로 명명하고 배정 가능
#     환자 수를 명시적으로 산출한다.
#
#  산출: results/16_lihc_stage_grade.csv
#        results/16_lihc_etiology.csv
#        results/16_lihc_serologynegative_cox.csv
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
hdr("16 · TCGA-LIHC 병기 · 조직등급 · 병인 층화")

## 1. 발현
c5 <- fetch_gene(XH$tcga, XDS$lihc_hiseq, "CYB5R3")
if (is.null(c5)) stop("CYB5R3 발현 조회 실패")
v <- tovec(c5); names(v) <- nb(names(v))
E <- data.frame(sample = names(v), CYB5R3 = as.numeric(v))
E <- E[typ(E$sample) %in% TUMOR_CODES, ]
E$patient <- pat(E$sample)
E <- E %>% dplyr::group_by(patient) %>% dplyr::slice(1) %>% dplyr::ungroup()

## 2. 임상
xq <- UCSCXenaTools::XenaGenerate() %>%
      UCSCXenaTools::XenaFilter(filterDatasets = XDS$lihc_clin)
xd <- UCSCXenaTools::XenaQuery(xq) %>%
      UCSCXenaTools::XenaDownload(destdir = CACHE, trans_slash = TRUE, force = FALSE)
cl <- UCSCXenaTools::XenaPrepare(xd)
if (is.list(cl) && !is.data.frame(cl)) cl <- cl[[1]]
cl <- as.data.frame(cl); names(cl)[1] <- "sample"; cl$patient <- pat(cl$sample)
pick <- function(p) { h <- grep(p, names(cl), ignore.case = TRUE, value = TRUE); if (length(h)) h[1] else NA }
c_stage <- pick("pathologic_stage|ajcc.*stage")
c_grade <- pick("neoplasm_histologic_grade|histologic.*grade")
## 바이러스 혈청검사 관련 열을 전부 찾는다 (표면항원/항체 등)
c_sero  <- grep("viral_hepatitis_serolog|hbv|hcv|hepatitis", names(cl),
                ignore.case = TRUE, value = TRUE)
message("  [cols] stage='", c_stage, "' grade='", c_grade, "'")
message("  [cols] serology: ", paste(c_sero, collapse = " | "))

cl2 <- cl %>% dplyr::group_by(patient) %>% dplyr::slice(1) %>% dplyr::ungroup()
D <- dplyr::inner_join(E[, c("patient","CYB5R3")], as.data.frame(cl2), by = "patient")

## 3. 병기·조직등급 경향성 (Jonckheere–Terpstra, 없으면 Kruskal–Wallis)
jt <- function(x, g) {
  g <- factor(g, levels = sort(unique(g)))
  if (requireNamespace("clinfun", quietly = TRUE))
    return(list(test = "Jonckheere-Terpstra",
                p = clinfun::jonckheere.test(x, as.numeric(g), alternative = "two.sided")$p.value))
  list(test = "Kruskal-Wallis (clinfun not installed)", p = kruskal.test(x, g)$p.value)
}
rows <- list()
for (nm in c(stage = c_stage, grade = c_grade)) {
  if (is.na(nm)) next
  lab <- names(which(c(stage = c_stage, grade = c_grade) == nm))[1]
  lv <- as.character(D[[nm]]); lv[!nzchar(trimws(lv))] <- NA
  lv <- if (lab == "stage")
    dplyr::case_when(grepl("Stage IV", lv) ~ "IV", grepl("Stage III", lv) ~ "III",
                     grepl("Stage II", lv) ~ "II", grepl("Stage I", lv) ~ "I", TRUE ~ NA_character_)
  else sub("^\\s*", "", lv)
  ok <- !is.na(lv) & !is.na(D$CYB5R3)
  if (sum(ok) < 30) next
  tt <- jt(D$CYB5R3[ok], lv[ok])
  tab <- table(lv[ok])
  med <- tapply(D$CYB5R3[ok], lv[ok], median)
  rows[[lab]] <- data.frame(
    stratifier = lab, test = tt$test, p = tt$p,
    levels = paste(sprintf("%s (n=%d, median %.2f)", names(tab), as.integer(tab),
                           med[names(tab)]), collapse = "; "))
}
SG <- dplyr::bind_rows(rows); print(SG); w_res(SG, "16_lihc_stage_grade.csv")

## 4. 병인 — TCGA-LIHC 의 'viral_hepatitis_serology' 는 **양성으로 나온 검사 항목을**
##    나열하는 열이다(예: 'Hepatitis B Surface Antigen'). 따라서
##      · 값이 있으면      = 바이러스 혈청검사 양성
##      · 값이 비어 있으면 = '기록 없음' 이며, 음성인지 미검사인지 구분 불가
##    리뷰어 3의 지적대로 후자는 비바이러스/대사성 병인과 동의어가 아니므로
##    'no positive viral serology recorded' 로만 명명한다.
D$sero <- NA_character_
if (length(c_sero)) {
  vals <- apply(sapply(c_sero, function(k) as.character(D[[k]])), 1, function(r) {
    r <- trimws(r[!is.na(r)])
    r <- r[nzchar(r) & !grepl("^\\[Not|^\\[Unknown|^NA$|^-$", r, ignore.case = TRUE)]
    paste(r, collapse = "; ")
  })
  message("  [diag] serology 원값 상위: ",
          paste(utils::head(sprintf("'%s'(%d)", names(sort(table(vals), decreasing = TRUE)),
                                    as.integer(sort(table(vals), decreasing = TRUE))), 6),
                collapse = " "))
  D$sero <- ifelse(nzchar(vals), "Viral serology positive",
                                 "No positive viral serology recorded")
}
tb <- table(D$sero, useNA = "ifany")
message(sprintf("  병인 분류: %s", paste(names(tb), as.integer(tb), collapse = " / ")))
ET <- data.frame(
  n_total = nrow(D),
  n_viral_positive = sum(D$sero == "Viral serology positive", na.rm = TRUE),
  n_no_positive_recorded = sum(D$sero == "No positive viral serology recorded", na.rm = TRUE))
if (ET$n_viral_positive >= 10 && ET$n_no_positive_recorded >= 10) {
  w <- wilcox.test(CYB5R3 ~ sero, data = D[!is.na(D$sero), ])
  ET$wilcox_p <- signif(w$p.value, 3)
  ET$median_viral_positive <- median(D$CYB5R3[D$sero == "Viral serology positive"], na.rm = TRUE)
  ET$median_no_positive <- median(D$CYB5R3[D$sero == "No positive viral serology recorded"], na.rm = TRUE)
}
ET$label_note <- paste("The TCGA-LIHC 'viral_hepatitis_serology' field lists serologies that were POSITIVE.",
  "An empty field means no positive serology was recorded and cannot be distinguished from 'not tested'.",
  "It is therefore NOT equivalent to non-viral or metabolic aetiology, and no metabolic subset can be defined from these data.")
print(ET); w_res(ET, "16_lihc_etiology.csv")

## 5. serology-negative 부분집합에서의 탐색적 생존분석 (음성 결과를 본문에 보고)
try({
  sv <- UCSCXenaTools::fetch_dense_values(XH$pan, XDS$pan_surv, c("OS","OS.time"),
                                          use_probeMap = FALSE)
  sv <- as.data.frame(t(sv)); sv$patient <- pat(rownames(sv))
  sv$OS <- suppressWarnings(as.integer(as.character(sv$OS)))
  sv$OS.time <- suppressWarnings(as.numeric(as.character(sv$OS.time)))
  d <- dplyr::inner_join(D[, c("patient","CYB5R3","sero")], sv, by = "patient")
  d <- d[!is.na(d$OS) & !is.na(d$OS.time) & d$OS.time > 0, ]
  out <- list()
  d <- d[!duplicated(d$patient), ]
  for (grp in c("No positive viral serology recorded","Viral serology positive")) {
    s <- d[which(d$sero == grp), ]
    if (nrow(s) < 10 || sum(s$OS == 1) < 3) next
    s$z <- as.numeric(scale(s$CYB5R3))
    fit <- summary(coxph(Surv(OS.time, OS) ~ z, data = s))
    out[[grp]] <- data.frame(subset = grp, n = nrow(s), events = sum(s$OS == 1),
                             HR_perSD = fit$conf.int[1,"exp(coef)"],
                             CI_low = fit$conf.int[1,"lower .95"],
                             CI_high = fit$conf.int[1,"upper .95"],
                             p = fit$coefficients[1,"Pr(>|z|)"],
                             note = "EXPLORATORY. 'No positive viral serology recorded' is NOT a metabolic subset.")
  }
  SN <- dplyr::bind_rows(out)
  if (nrow(SN)) { print(SN); w_res(SN, "16_lihc_serologynegative_cox.csv") }
}, silent = FALSE)

save_session("16_lihc_stratified")
message("\n[16] 완료.")
