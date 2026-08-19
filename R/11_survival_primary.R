# =====================================================================
#  11_survival_primary.R   —  PRIMARY PROGNOSTIC ANALYSIS
#  목적 (심사 대응):
#   · R1-4 "The continuous multivariable Cox model should be designated as the
#     primary prognostic analysis" + 사건수·제외·추적기간·결측처리·병기부호화·
#     PH 진단·비선형성 보고
#   · R3 "harden the survival modelling": 환자 흐름(371→339), 사망수,
#     PH 가정 검정, C-index 및 병기·연령·성별 기저모형 대비 ΔC-index,
#     AFP·혈관침습·간경변 공변량 추가
#   · R1-3 "include both genes in the same multivariable survival model, and
#     test whether a combined or interaction model adds information"
#   · R1-4 / R3 "GSE14520 … whether adjusted for clinical covariates"
#
#  산출: results/11_flow.csv                 환자 흐름 + 사건수 + 추적기간
#        results/11_cox_primary.csv          주분석(연속형 per-SD)
#        results/11_cox_ph.csv               PH 가정 검정
#        results/11_cox_nonlinearity.csv     스플라인 대 선형 LRT
#        results/11_cindex.csv               C-index 및 ΔC (부트스트랩 CI)
#        results/11_cox_joint.csv            결합·상호작용 모형 + LRT
#        results/11_cox_extended.csv         AFP/혈관침습/섬유화 보정 모형
#        results/11_missing_compare.csv      포함 vs 제외 환자 비교
#        results/11_cox_imputed.csv          다중대치 민감도분석
#        results/11_gse14520_adjusted.csv    외부 코호트 보정 모형
#        figures/11_ph_schoenfeld.png, figures/11_forest_primary.png
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
hdr("11 · 주 예후분석 (TCGA-LIHC) + 외부검증")

flow <- list()
add_flow <- function(step, n, note = "") {
  flow[[length(flow)+1]] <<- data.frame(step = step, n = n, note = note)
  message(sprintf("  [flow] %-52s n = %s", step, n))
}

## =====================================================================
## 1. 발현
## =====================================================================
c5 <- fetch_gene(XH$tcga, XDS$lihc_hiseq, "CYB5R3")
if (is.null(c5)) stop("CYB5R3 발현 조회 실패 (HiSeqV2)")
v5 <- tovec(c5); names(v5) <- nb(names(v5))
C <- data.frame(sample = names(v5), CYB5R3 = as.numeric(v5))
add_flow("TCGA-LIHC HiSeqV2 samples with CYB5R3", nrow(C))
C <- C[typ(C$sample) %in% TUMOR_CODES, ]
C$patient <- pat(C$sample)
C <- C %>% dplyr::group_by(patient) %>% dplyr::slice(1) %>% dplyr::ungroup()
add_flow("Primary tumors, one sample per patient (CYB5R3)", nrow(C))

## MTARC1 은 legacy probeMap 에 없을 수 있어 TOIL 허브를 사용 (원 투고와 동일)
mt <- fetch_gene(XH$toil, XDS$toil_tpm, "MTARC1")
if (is.null(mt)) stop("MTARC1 발현 조회 실패 (TOIL)")
vm <- tovec(mt); names(vm) <- nb(names(vm))
M <- data.frame(sample = names(vm), MTARC1 = as.numeric(vm))
M <- M[grepl("^TCGA", M$sample) & typ(M$sample) %in% TUMOR_CODES, ]
M$patient <- pat(M$sample)
M <- M %>% dplyr::group_by(patient) %>% dplyr::slice(1) %>% dplyr::ungroup()
E <- dplyr::inner_join(C[, c("sample","patient","CYB5R3")],
                       M[, c("patient","MTARC1")], by = "patient")
add_flow("Merged CYB5R3 + MTARC1 (patient level)", nrow(E),
         sprintf("MTARC1 alias used: %s", attr(mt, "alias_used")))

## =====================================================================
## 2. 임상 (범주값 보존을 위해 clinicalMatrix 를 표로 내려받음)
## =====================================================================
xq <- UCSCXenaTools::XenaGenerate() %>%
      UCSCXenaTools::XenaFilter(filterDatasets = XDS$lihc_clin)
xd <- UCSCXenaTools::XenaQuery(xq) %>%
      UCSCXenaTools::XenaDownload(destdir = CACHE, trans_slash = TRUE, force = FALSE)
cl <- UCSCXenaTools::XenaPrepare(xd)
if (is.list(cl) && !is.data.frame(cl)) cl <- cl[[1]]
cl <- as.data.frame(cl); names(cl)[1] <- "sample"; cl$patient <- pat(cl$sample)

pick <- function(p) { h <- grep(p, names(cl), ignore.case = TRUE, value = TRUE); if (length(h)) h[1] else NA }
col_age  <- pick("^age_at_initial|^age$")
col_sex  <- pick("^gender$|^sex$")
col_stg  <- pick("pathologic_stage|^stage$|ajcc.*stage")
col_grd  <- pick("neoplasm_histologic_grade|histologic.*grade")
## AFP: 'fetoprotein_outcome_value' 가 실제 측정치이고
##      '..._lower_limit' / '..._upper_limit' 는 검사 한계값이므로 value 를 우선한다.
col_afp  <- pick("fetoprotein.*value|^afp.*value|^afp$")
if (is.na(col_afp)) col_afp <- pick("fetoprotein|^afp")
col_vasc <- pick("vascular_tumor_cell_type|vascular.*invasion")
col_fib  <- pick("fibrosis_ishak|ishak")
col_cp   <- pick("child_pugh")
message(sprintf("  [cols] age='%s' sex='%s' stage='%s' grade='%s'\n         AFP='%s' vascular='%s' Ishak='%s' ChildPugh='%s'",
                col_age, col_sex, col_stg, col_grd, col_afp, col_vasc, col_fib, col_cp))

cl$age   <- suppressWarnings(as.numeric(as.character(cl[[col_age]])))
cl$sex   <- as.character(cl[[col_sex]])
cl$stage_raw <- as.character(cl[[col_stg]])
## 병기 부호화 명시 (심사 R1-4): AJCC I/II = 'I-II', III(A/B/C)/IV(A/B) = 'III-IV'
cl$stage_bin <- dplyr::case_when(
  stringr::str_detect(cl$stage_raw, "Stage III|Stage IV") ~ "III-IV",
  stringr::str_detect(cl$stage_raw, "Stage I{1,2}$|Stage I |Stage II") ~ "I-II",
  stringr::str_detect(cl$stage_raw, "^Stage I$")          ~ "I-II",
  TRUE ~ NA_character_)
cl$grade <- if (!is.na(col_grd)) as.character(cl[[col_grd]]) else NA_character_
cl$afp   <- if (!is.na(col_afp)) suppressWarnings(as.numeric(as.character(cl[[col_afp]]))) else NA_real_
cl$vasc  <- if (!is.na(col_vasc)) as.character(cl[[col_vasc]]) else NA_character_
cl$ishak <- if (!is.na(col_fib)) as.character(cl[[col_fib]]) else NA_character_
## 간경변 정의: Ishak 5–6 = 확립된 간경변
cl$cirrhosis <- ifelse(is.na(cl$ishak), NA,
  ifelse(grepl("^5|^6", trimws(cl$ishak)), "Cirrhosis (Ishak 5-6)", "No cirrhosis (Ishak 0-4)"))
cl$vasc_bin <- ifelse(is.na(cl$vasc), NA,
  ifelse(grepl("None", cl$vasc, ignore.case = TRUE), "None", "Micro/Macro"))
cl <- cl %>% dplyr::group_by(patient) %>% dplyr::slice(1) %>% dplyr::ungroup()

## =====================================================================
## 3. 생존 (TCGA-CDR)
## =====================================================================
sv <- UCSCXenaTools::fetch_dense_values(
  XH$pan, XDS$pan_surv, c("OS","OS.time","PFI","PFI.time"), use_probeMap = FALSE)
sv <- as.data.frame(t(sv)); sv$patient <- pat(rownames(sv))
for (k in c("OS","OS.time","PFI","PFI.time"))
  sv[[k]] <- suppressWarnings(as.numeric(as.character(sv[[k]])))
sv <- sv %>% dplyr::group_by(patient) %>% dplyr::slice(1) %>% dplyr::ungroup()

dat <- E %>%
  dplyr::left_join(cl[, c("patient","age","sex","stage_raw","stage_bin","grade",
                          "afp","vasc_bin","cirrhosis")], by = "patient") %>%
  dplyr::left_join(sv[, c("patient","OS","OS.time","PFI","PFI.time")], by = "patient")
add_flow("With clinical + survival annotation", nrow(dat))
dat_os <- dat %>% dplyr::filter(!is.na(OS), !is.na(OS.time), OS.time > 0)
add_flow("With non-missing OS time and status", nrow(dat_os))
cc <- dat_os %>% dplyr::filter(!is.na(age), !is.na(sex), !is.na(stage_bin))
add_flow("PRIMARY ANALYSIS SET (complete age, sex, stage)", nrow(cc),
         sprintf("deaths = %d", sum(cc$OS == 1)))
mfu <- median_followup(cc$OS.time, cc$OS)
add_flow("Median follow-up, days (reverse Kaplan-Meier)", round(mfu, 1),
         sprintf("= %.1f months", mfu / 30.44))

cc$z5 <- as.numeric(scale(cc$CYB5R3))
cc$zm <- as.numeric(scale(cc$MTARC1))
cc$sex <- factor(cc$sex)
cc$stage_bin <- factor(cc$stage_bin, levels = c("I-II","III-IV"))
w_res(dplyr::bind_rows(flow), "11_flow.csv")

## =====================================================================
## 4. 주분석 — 연속형(per-SD) 다변량 Cox
## =====================================================================
hdr("4. PRIMARY: continuous per-SD multivariable Cox (OS)")
tidy_cox <- function(fit, label, term) {
  s <- summary(fit); i <- which(rownames(s$coefficients) == term)
  data.frame(model = label, term = term, n = s$n, events = s$nevent,
             HR = s$conf.int[i,"exp(coef)"],
             CI_low = s$conf.int[i,"lower .95"], CI_high = s$conf.int[i,"upper .95"],
             p = s$coefficients[i,"Pr(>|z|)"])
}
prim <- list()
f5 <- coxph(Surv(OS.time, OS) ~ z5 + age + sex + stage_bin, data = cc)
fm <- coxph(Surv(OS.time, OS) ~ zm + age + sex + stage_bin, data = cc)
prim[["CYB5R3"]] <- tidy_cox(f5, "PRIMARY: OS ~ CYB5R3(per SD) + age + sex + stage", "z5")
prim[["MTARC1"]] <- tidy_cox(fm, "PRIMARY: OS ~ MTARC1(per SD) + age + sex + stage", "zm")
## 사전지정 민감도 — PFI 종점, 그리고 중앙값 이분(참고용, 주분석 아님)
pf <- dat_os %>% dplyr::filter(!is.na(PFI), !is.na(PFI.time), PFI.time > 0,
                               !is.na(age), !is.na(sex), !is.na(stage_bin))
pf$z5 <- as.numeric(scale(pf$CYB5R3)); pf$zm <- as.numeric(scale(pf$MTARC1))
pf$sex <- factor(pf$sex); pf$stage_bin <- factor(pf$stage_bin, c("I-II","III-IV"))
prim[["CYB5R3_PFI"]] <- tidy_cox(coxph(Surv(PFI.time, PFI) ~ z5 + age + sex + stage_bin, data = pf),
                                 "SENSITIVITY: PFI ~ CYB5R3(per SD) + covariates", "z5")
prim[["MTARC1_PFI"]] <- tidy_cox(coxph(Surv(PFI.time, PFI) ~ zm + age + sex + stage_bin, data = pf),
                                 "SENSITIVITY: PFI ~ MTARC1(per SD) + covariates", "zm")
cc$hi5 <- factor(ifelse(cc$CYB5R3 >= median(cc$CYB5R3), "High","Low"), c("Low","High"))
prim[["CYB5R3_median"]] <- tidy_cox(coxph(Surv(OS.time, OS) ~ hi5 + age + sex + stage_bin, data = cc),
                                    "SECONDARY (not primary): median split", "hi5High")
res_prim <- dplyr::bind_rows(prim); print(res_prim); w_res(res_prim, "11_cox_primary.csv")

## =====================================================================
## 5. 비례위험(PH) 가정
## =====================================================================
hdr("5. Proportional-hazards diagnostics")
ph5 <- cox.zph(f5); phm <- cox.zph(fm)
ph_df <- rbind(
  data.frame(model = "CYB5R3", term = rownames(ph5$table), chisq = ph5$table[,"chisq"],
             df = ph5$table[,"df"], p = ph5$table[,"p"]),
  data.frame(model = "MTARC1", term = rownames(phm$table), chisq = phm$table[,"chisq"],
             df = phm$table[,"df"], p = phm$table[,"p"]))
print(ph_df); w_res(ph_df, "11_cox_ph.csv")
png(file.path(FIG, "11_ph_schoenfeld.png"), width = 1600, height = 1200, res = 150)
par(mfrow = c(2,2)); try(plot(ph5)); dev.off()

## =====================================================================
## 6. 비선형 발현효과 (심사 R1-4 "possible nonlinear expression effects")
## =====================================================================
hdr("6. Nonlinearity of the expression effect")
nl <- tryCatch({
  fit_lin <- coxph(Surv(OS.time, OS) ~ z5 + age + sex + stage_bin, data = cc)
  fit_spl <- coxph(Surv(OS.time, OS) ~ pspline(z5, df = 3) + age + sex + stage_bin, data = cc)
  lr <- 2 * (fit_spl$loglik[2] - fit_lin$loglik[2])
  ddf <- abs(diff(c(fit_lin$df %||% length(coef(fit_lin)),
                    sum(fit_spl$df))))
  data.frame(gene = "CYB5R3", test = "linear vs penalised spline (df=3)",
             LR_chisq = lr, df = ddf, p = pchisq(lr, max(ddf, 1), lower.tail = FALSE),
             AIC_linear = AIC(fit_lin), AIC_spline = AIC(fit_spl))
}, error = function(e) data.frame(gene = "CYB5R3", test = "spline failed",
                                  LR_chisq = NA, df = NA, p = NA,
                                  AIC_linear = NA, AIC_spline = NA))
print(nl); w_res(nl, "11_cox_nonlinearity.csv")

## =====================================================================
## 7. C-index 와 ΔC-index (심사 R3: "does CYB5R3 add anything to staging?")
## =====================================================================
hdr("7. Discrimination: C-index and change over the clinical base model")
ci5 <- c_index_delta(cc, "OS.time", "OS", "age + sex + stage_bin",
                     "age + sex + stage_bin + z5", B = 1000)
cim <- c_index_delta(cc, "OS.time", "OS", "age + sex + stage_bin",
                     "age + sex + stage_bin + zm", B = 1000)
cij <- c_index_delta(cc, "OS.time", "OS", "age + sex + stage_bin",
                     "age + sex + stage_bin + z5 + zm", B = 1000)
cidx <- dplyr::bind_rows(cbind(added = "CYB5R3", ci5), cbind(added = "MTARC1", cim),
                         cbind(added = "CYB5R3 + MTARC1", cij))
print(cidx); w_res(cidx, "11_cindex.csv")

## =====================================================================
## 8. 결합 모형과 상호작용 (심사 R1-3)
## =====================================================================
hdr("8. Joint and interaction models for the two axis members")
f_base <- coxph(Surv(OS.time, OS) ~ age + sex + stage_bin, data = cc)
f_j    <- coxph(Surv(OS.time, OS) ~ z5 + zm + age + sex + stage_bin, data = cc)
f_i    <- coxph(Surv(OS.time, OS) ~ z5 * zm + age + sex + stage_bin, data = cc)
## 축 불균형(axis imbalance) 지표: CYB5R3 공급 − mARC1 수용, 표준화 차
cc$axis_gap <- cc$z5 - cc$zm
f_g    <- coxph(Surv(OS.time, OS) ~ axis_gap + age + sex + stage_bin, data = cc)
lrt <- function(a, b, lab) {
  d <- 2 * (b$loglik[2] - a$loglik[2]); df <- length(coef(b)) - length(coef(a))
  data.frame(comparison = lab, LR_chisq = d, df = df,
             p = pchisq(d, max(df,1), lower.tail = FALSE))
}
joint <- dplyr::bind_rows(
  tidy_cox(f_j, "JOINT: OS ~ CYB5R3 + MTARC1 + covariates", "z5"),
  tidy_cox(f_j, "JOINT: OS ~ CYB5R3 + MTARC1 + covariates", "zm"),
  tidy_cox(f_i, "INTERACTION: OS ~ CYB5R3 * MTARC1 + covariates", "z5:zm"),
  tidy_cox(f_g, "AXIS GAP: OS ~ (z_CYB5R3 - z_MTARC1) + covariates", "axis_gap"))
joint_lrt <- dplyr::bind_rows(
  lrt(f_base, f5, "base + CYB5R3 vs base"),
  lrt(f_base, f_j, "base + CYB5R3 + MTARC1 vs base"),
  lrt(f5,  f_j, "joint vs CYB5R3 alone (does MTARC1 add?)"),
  lrt(f_j, f_i, "interaction vs joint (is the axis multiplicative?)"))
print(joint); print(joint_lrt)
w_res(joint, "11_cox_joint.csv")
w_res(joint_lrt, "11_cox_joint_lrt.csv")
## CYB5R3–MTARC1 표본수준 상관 (심사 R1-3 "examine their sample-level relationship")
sl <- spearman_ci(cc$CYB5R3, cc$MTARC1)
sl$comparison <- "CYB5R3 vs MTARC1, TCGA-LIHC primary tumors (Spearman)"
w_res(sl, "11_axis_sample_level_correlation.csv"); print(sl)

## =====================================================================
## 9. 확장 공변량 (AFP · 혈관침습 · 간경변)
## =====================================================================
hdr("9. Extended-covariate models (AFP, vascular invasion, cirrhosis)")
ext <- list()
addm <- function(extra, label) {
  d <- cc
  for (v in extra) d <- d[!is.na(d[[v]]), ]
  if (nrow(d) < 60 || sum(d$OS == 1) < 20) {
    message("  [skip] ", label, " — n=", nrow(d), " events=", sum(d$OS == 1)); return(NULL)
  }
  for (v in extra) if (is.character(d[[v]])) d[[v]] <- factor(d[[v]])
  d$z5 <- as.numeric(scale(d$CYB5R3))
  f <- as.formula(paste("Surv(OS.time, OS) ~ z5 + age + sex + stage_bin +",
                        paste(extra, collapse = " + ")))
  tidy_cox(coxph(f, data = d), label, "z5")
}
push <- function(x) if (!is.null(x)) ext[[length(ext)+1]] <<- x
if (!all(is.na(cc$afp))) {
  cc$log_afp <- log10(pmax(cc$afp, 1))
  push(tryCatch(addm("log_afp", "EXTENDED: + log10(AFP)"), error = function(e) NULL))
} else message("  [skip] AFP: 임상표에 값이 없음")
push(tryCatch(addm("vasc_bin",  "EXTENDED: + vascular invasion"), error = function(e) NULL))
push(tryCatch(addm("cirrhosis", "EXTENDED: + cirrhosis (Ishak 5-6)"), error = function(e) NULL))
push(tryCatch(addm(c("vasc_bin","cirrhosis"), "EXTENDED: + vascular invasion + cirrhosis"),
              error = function(e) NULL))
ext <- dplyr::bind_rows(ext)
if (nrow(ext)) { print(ext); w_res(ext, "11_cox_extended.csv") } else
  message("  확장 공변량이 임상표에 없거나 결측이 많아 모형을 적합하지 못했습니다 (원고에 그대로 기술).")
## 어떤 확장 공변량이 실제로 존재/가용했는지 기록
avail <- data.frame(
  covariate = c("AFP","vascular invasion","Ishak fibrosis / cirrhosis","histologic grade","Child-Pugh"),
  column    = c(col_afp, col_vasc, col_fib, col_grd, col_cp),
  n_nonmissing = c(sum(!is.na(cc$afp)), sum(!is.na(cc$vasc_bin)),
                   sum(!is.na(cc$cirrhosis)), sum(!is.na(cc$grade)),
                   if (!is.na(col_cp)) sum(!is.na(cl[[col_cp]])) else 0))
w_res(avail, "11_extended_covariate_availability.csv"); print(avail)

## =====================================================================
## 10. 결측 처리 — 포함 vs 제외 비교, 다중대치 민감도
## =====================================================================
hdr("10. Missing data: included vs excluded, and multiple imputation")
dat_os$included <- dat_os$patient %in% cc$patient
cmp <- dat_os %>% dplyr::group_by(included) %>% dplyr::summarise(
  n = dplyr::n(), deaths = sum(OS == 1),
  median_CYB5R3 = median(CYB5R3, na.rm = TRUE),
  median_MTARC1 = median(MTARC1, na.rm = TRUE),
  median_age = median(age, na.rm = TRUE),
  pct_male = 100 * mean(sex == "MALE", na.rm = TRUE),
  pct_stage_III_IV = 100 * mean(stage_bin == "III-IV", na.rm = TRUE),
  median_OS_days = median(OS.time, na.rm = TRUE), .groups = "drop")
pw <- suppressWarnings(wilcox.test(CYB5R3 ~ included, data = dat_os)$p.value)
cmp$wilcox_p_CYB5R3_included_vs_excluded <- pw
print(cmp); w_res(cmp, "11_missing_compare.csv")

if (has_pkg("mice")) {
  mi <- tryCatch({
    d <- dat_os[, c("CYB5R3","MTARC1","age","sex","stage_bin","OS","OS.time")]
    d$sex <- factor(d$sex); d$stage_bin <- factor(d$stage_bin, c("I-II","III-IV"))
    ## 생존자료 대치의 표준: Nelson–Aalen 누적위험을 보조변수로 포함
    bh <- survival::basehaz(coxph(Surv(OS.time, OS) ~ 1, data = d))
    d$na_cumhaz <- bh$hazard[pmax(findInterval(d$OS.time, bh$time), 1)]
    imp <- mice::mice(d, m = 20, seed = 20260817, printFlag = FALSE)
    fitm <- with(imp, coxph(Surv(OS.time, OS) ~ scale(CYB5R3)[,1] + age + sex + stage_bin))
    po <- summary(mice::pool(fitm), conf.int = TRUE)
    po <- po[1, ]
    data.frame(model = "MULTIPLE IMPUTATION (m=20): OS ~ CYB5R3(per SD) + covariates",
               n = nrow(d), events = sum(d$OS == 1),
               HR = exp(po$estimate), CI_low = exp(po$`2.5 %`), CI_high = exp(po$`97.5 %`),
               p = po$p.value)
  }, error = function(e) { message("  [mice 실패] ", conditionMessage(e)); NULL })
  if (!is.null(mi)) { print(mi); w_res(mi, "11_cox_imputed.csv") }
} else message("  mice 패키지가 없어 다중대치 민감도분석을 건너뜁니다.")

## 포레스트 그림
try({
  fp <- res_prim[grepl("^PRIMARY", res_prim$model), ]
  fp$lab <- ifelse(fp$term == "z5", "CYB5R3", "MTARC1")
  g <- ggplot2::ggplot(fp, ggplot2::aes(HR, lab)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey55") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = CI_low, xmax = CI_high), height = .12) +
    ggplot2::geom_point(size = 3, colour = "#1F3A5F") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("HR %.2f (%.2f-%.2f), p=%.3g", HR, CI_low, CI_high, p)),
                       vjust = -1.1, size = 3) +
    ggplot2::scale_x_log10() + ggplot2::theme_minimal() +
    ggplot2::labs(x = "Adjusted HR for overall survival per SD (age, sex, AJCC stage)",
                  y = NULL, title = sprintf("TCGA-LIHC primary analysis (n=%d, %d deaths)",
                                            nrow(cc), sum(cc$OS == 1)))
  ggplot2::ggsave(file.path(FIG, "11_forest_primary.png"), g, width = 7, height = 3, dpi = 300)
}, silent = TRUE)

## =====================================================================
## 11. 외부 코호트 GSE14520 — 전처리·프로브 선택 명시 + 공변량 보정
## =====================================================================
hdr("11. External cohort GSE14520 — documented preprocessing and adjusted models")
try({
  library(GEOquery); library(Biobase)
  gl <- getGEO("GSE14520", GSEMatrix = TRUE, destdir = CACHE)
  message("  ExpressionSets: ", length(gl))
  em <- list(); probe_log <- list()
  for (i in seq_along(gl)) {
    es <- gl[[i]]; ex <- exprs(es); fd <- fData(es)
    scol <- grep("Gene.?Symbol|GENE_SYMBOL|Symbol", names(fd), ignore.case = TRUE, value = TRUE)[1]
    if (is.na(scol)) next
    syms <- toupper(as.character(fd[[scol]]))
    plat <- annotation(es)
    row <- data.frame(geo = colnames(ex))
    for (g in AXIS) {
      idx <- which(syms %in% toupper(ALIAS[[g]]))
      if (!length(idx)) next
      ## 프로브 선택 규칙 (사전지정): 해당 심볼에 매핑된 모든 프로브의 평균.
      ## 프로브 ID 와 개수를 기록하여 재현 가능하게 한다.
      row[[g]] <- as.numeric(colMeans(ex[idx, , drop = FALSE], na.rm = TRUE))
      probe_log[[length(probe_log)+1]] <- data.frame(
        platform = plat, gene = g, n_probes = length(idx),
        probes = paste(rownames(ex)[idx], collapse = ";"),
        rule = "mean of all probes mapped to the symbol (pre-specified)",
        preprocessing = "GEO series matrix as deposited (RMA-normalised, log2); no further transformation")
    }
    if (ncol(row) > 1) em[[length(em)+1]] <- row
  }
  em <- dplyr::bind_rows(em)
  w_res(dplyr::bind_rows(probe_log), "11_gse14520_probe_selection.csv")

  ## 임상 부속표
  cl2 <- NULL
  for (u in c("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/suppl/GSE14520_Extra_Supplement.txt.gz",
              "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/suppl/GSE14520_Extra_Supplement.txt")) {
    f <- file.path(CACHE, basename(u))
    ok <- file.exists(f) || tryCatch({ download.file(u, f, mode = "wb", quiet = TRUE); TRUE },
                                     error = function(e) FALSE)
    if (ok) { cl2 <- tryCatch(data.table::fread(f), error = function(e) NULL); if (!is.null(cl2)) break }
  }
  stopifnot(!is.null(cl2))
  gg  <- function(p) grep(p, names(cl2), ignore.case = TRUE, value = TRUE)[1]
  key <- data.frame(
    geo    = as.character(cl2[[gg("Affy.?GSM")]]),
    time   = suppressWarnings(as.numeric(cl2[[gg("Survival.*month")]])),
    status = suppressWarnings(as.integer(cl2[[gg("Survival.*status")]])),
    sexv   = as.character(cl2[[gg("^Gender")]]),
    agev   = suppressWarnings(as.numeric(cl2[[gg("^Age")]])),
    cirrh  = as.character(cl2[[gg("Cirrhosis")]]),
    tnm    = as.character(cl2[[gg("TNM")]]),
    bclc   = as.character(cl2[[gg("BCLC")]]),
    afp    = as.character(cl2[[gg("AFP")]]),
    size   = as.character(cl2[[gg("Tumor Size")]]),
    multin = as.character(cl2[[gg("Multinodular")]]),
    tissue = as.character(cl2[[gg("Tissue Type")]]))
  m <- merge(em, key, by = "geo")
  ## 종점 정의 (심사 R1-4): 전체생존, 개월 단위, 사망=1
  m <- m[!is.na(m$time) & m$time > 0 & !is.na(m$status), ]
  message("  merged n = ", nrow(m), " ; deaths = ", sum(m$status == 1))
  ## 공변량 원값 분포를 먼저 찍는다 (부호화 실패를 눈으로 확인하기 위해)
  for (v in c("sexv","cirrh","tnm","bclc","afp","size","multin","tissue"))
    message(sprintf("  [diag] %-7s : %s", v,
      paste(utils::head(sprintf("%s(%d)", names(table(m[[v]], useNA = "ifany")),
                                as.integer(table(m[[v]], useNA = "ifany"))), 8), collapse = " ")))
  tn <- toupper(trimws(as.character(m$tnm)))
  m$tnm_bin <- ifelse(grepl("^(III|IV|3|4)", tn), "III-IV",
               ifelse(grepl("^(I|II|1|2)",  tn), "I-II", NA))
  m$afp_bin  <- ifelse(is.na(m$afp),  NA, ifelse(grepl(">", m$afp),  "high",  "low"))
  m$size_bin <- ifelse(is.na(m$size), NA, ifelse(grepl(">", m$size), "large", "small"))
  message("  [diag] tnm_bin: ", paste(names(table(m$tnm_bin, useNA = "ifany")),
                                      as.integer(table(m$tnm_bin, useNA = "ifany")),
                                      collapse = " / "))

  ## 사용 가능한 공변량만 골라 모형을 세운다 (수준이 1개면 자동 제외)
  build <- function(d, vars) {
    keep <- character(0)
    for (v in vars) {
      if (!v %in% names(d)) next
      if (is.numeric(d[[v]])) { if (length(unique(na.omit(d[[v]]))) > 2) keep <- c(keep, v) }
      else { f <- safe_factor(d[[v]]); if (!is.null(f)) { d[[v]] <- f; keep <- c(keep, v) } }
    }
    dropped <- setdiff(vars, keep)
    if (length(dropped))
      message("  [covariate dropped] ", paste(dropped, collapse = ", "),
              " — absent from the phenotype table or single-level after filtering")
    list(d = d, keep = keep, dropped = dropped)
  }
  fit_one <- function(d, vars, gene, label) {
    d <- d[!is.na(d[[gene]]), ]
    for (v in vars) if (v %in% names(d)) d <- d[!is.na(d[[v]]), ]
    b <- build(d, vars); d <- b$d
    if (nrow(d) < 60 || sum(d$status == 1) < 20) {
      message("  [skip] ", label, " — n=", nrow(d), " events=", sum(d$status == 1)); return(NULL) }
    d$z <- as.numeric(scale(d[[gene]]))
    rhs <- paste(c("z", b$keep), collapse = " + ")
    tryCatch(tidy_cox(coxph(as.formula(paste("Surv(time, status) ~", rhs)), data = d),
                      paste0(label, "  [terms: ", rhs, "]"), "z"),
             error = function(e) { message("  [fail] ", label, ": ", conditionMessage(e)); NULL })
  }
  ext14 <- list()
  for (g in AXIS) {
    if (!g %in% names(m)) next
    ext14[[paste0(g,"_unadj")]] <- fit_one(m, character(0), g,
      sprintf("GSE14520 UNADJUSTED: OS ~ %s (per SD)", g))
    ext14[[paste0(g,"_adj")]]   <- fit_one(m, c("agev","sexv","tnm_bin"), g,
      sprintf("GSE14520 ADJUSTED (age, sex, TNM): OS ~ %s (per SD)", g))
    ## 라벨은 '요청한' 공변량이 아니라 '실제로 적합된' 항을 반영해야 한다.
    ## (이전 판은 AFP 를 라벨에 넣었지만 GSE14520 표현형표에 없어 모형에서 빠졌다.)
    ext14[[paste0(g,"_full")]]  <- fit_one(m, c("agev","sexv","tnm_bin","cirrh","afp_bin","multin"), g,
      sprintf("GSE14520 EXTENDED ADJUSTMENT: OS ~ %s (per SD)", g))
  }
  ## 두 유전자 결합 모형
  try({
    d4 <- m[!is.na(m$CYB5R3) & !is.na(m$MTARC1) & !is.na(m$agev) & !is.na(m$sexv) & !is.na(m$tnm_bin), ]
    b <- build(d4, c("agev","sexv","tnm_bin")); d4 <- b$d
    if (nrow(d4) > 80) {
      d4$z5 <- as.numeric(scale(d4$CYB5R3)); d4$zm <- as.numeric(scale(d4$MTARC1))
      rhs <- paste(c("z5","zm", b$keep), collapse = " + ")
      fj <- coxph(as.formula(paste("Surv(time, status) ~", rhs)), data = d4)
      ext14[["joint5"]] <- tidy_cox(fj, "GSE14520 JOINT (adjusted): CYB5R3 + MTARC1", "z5")
      ext14[["jointm"]] <- tidy_cox(fj, "GSE14520 JOINT (adjusted): CYB5R3 + MTARC1", "zm")
      d4$axis_gap <- d4$z5 - d4$zm
      ext14[["gap"]] <- tidy_cox(
        coxph(as.formula(paste("Surv(time, status) ~ axis_gap +", paste(b$keep, collapse = " + "))), data = d4),
        "GSE14520 AXIS GAP (adjusted): (z_CYB5R3 - z_MTARC1)", "axis_gap")
    }
  }, silent = FALSE)
  r14 <- dplyr::bind_rows(ext14); print(r14); w_res(r14, "11_gse14520_adjusted.csv")
}, silent = FALSE)

save_session("11_survival")
message("\n[11] 완료.")
