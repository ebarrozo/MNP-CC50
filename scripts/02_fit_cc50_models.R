# ==============================================================================
# mnp_CC50_analysis_v6_3_FINAL_BEWO.R
# Final NanoImpact CC50 analysis with BeWo as biologically motivated reference
#
# Primary model:
#   cc50 ~ cell_line + Treatment + Size + Timepoint + (1 | BioRep_ID)
#
# Secondary one-at-a-time interaction models:
#   cell_line * Treatment + Size + Timepoint + (1 | BioRep_ID)
#   cell_line * Timepoint + Size + Treatment + (1 | BioRep_ID)
#   cell_line * Size + Timepoint + Treatment + (1 | BioRep_ID)
#
# Sensitivities:
#   - Strict SSMD >= 1 additive LMM
#   - Additive ordinal CLMM:
#       <=0.2 < 0.2-20 < >20
#
# Cell-line presentation:
#   - BeWo-referenced adjusted contrasts (primary manuscript contrasts)
#   - All pairwise adjusted cell-line contrasts (supplementary)
#
# The model fit and omnibus tests do NOT depend on the reference category.
# ==============================================================================

# Run this script from the repository root.
project_dir <- normalizePath(getwd(), mustWork = TRUE)
if (!dir.exists(file.path(project_dir, "scripts")) ||
    !dir.exists(file.path(project_dir, "metadata"))) {
  stop("Run this script from the repository root (the directory containing scripts/ and metadata/).")
}


suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

if (!requireNamespace("ordinal", quietly = TRUE)) {
  stop("Package 'ordinal' is required. Install with install.packages('ordinal').")
}

# ------------------------------------------------------------------------------
# 0. Paths/settings
# ------------------------------------------------------------------------------

v5_dir <- file.path(project_dir, "results", "v5_2")
out_dir <- file.path(project_dir, "results", "v6_3")
plot_dir <- file.path(out_dir, "plots")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

primary_file <- file.path(v5_dir, "CC50_ReplicateValues.csv")
strict_file  <- file.path(v5_dir, "CC50_ReplicateValues_MEDIA_STRICT_SSMD_SENSITIVITY_ALL480.csv")
plate_qc_file <- file.path(v5_dir, "MNP_v5_2_PlateControl_QC.csv")

for (f in c(primary_file, strict_file, plate_qc_file)) {
  if (!file.exists(f)) stop("Missing required input: ", f)
}

LOW_CAP  <- 0.2
HIGH_CAP <- 20
REF_CELL <- "BeWo"
ALPHA <- 0.05

cell_levels <- c("BeWo", "HTR8", "JEG3", "THP-1", "Jurkat")
size_levels <- c("Small", "Large")
time_levels <- c("3", "6", "24", "48")
treatment_levels <- c("Mix", "PE", "PS", "PVC")

write_out <- function(df, filename) {
  readr::write_csv(df, file.path(out_dir, filename))
}

# ------------------------------------------------------------------------------
# 1. Standardize data
# ------------------------------------------------------------------------------

standardize_cc50 <- function(df) {
  if (!"cc50_censor" %in% names(df) && "CC50_Censor" %in% names(df)) {
    df <- df %>% rename(cc50_censor = CC50_Censor)
  }

  req <- c("cell_line","Size","Timepoint","Treatment","Replicate","cc50","cc50_censor")
  miss <- setdiff(req, names(df))
  if (length(miss) > 0) stop("Missing CC50 columns: ", paste(miss, collapse=", "))

  df %>%
    mutate(
      cell_line = recode(as.character(cell_line),
                         "THP1"="THP-1", "Jeg3"="JEG3",
                         .default=as.character(cell_line)),
      Size = recode(as.character(Size),
                    "small"="Small", "large"="Large",
                    .default=as.character(Size)),
      Timepoint = stringr::str_extract(as.character(Timepoint), "\\d+"),
      Treatment = as.character(Treatment),
      Replicate = as.character(Replicate),
      cc50 = suppressWarnings(as.numeric(cc50)),
      BioRep_ID = interaction(cell_line, Replicate, drop=TRUE),
      cell_line = factor(cell_line, levels=cell_levels),
      Treatment = factor(Treatment, levels=treatment_levels),
      Size = factor(Size, levels=size_levels),
      Timepoint = factor(Timepoint, levels=time_levels),
      Censor_Category = case_when(
        str_detect(toupper(as.character(cc50_censor)), "LOWER|<=|BELOW") ~ "Lower",
        str_detect(toupper(as.character(cc50_censor)), "UPPER|>|ABOVE") ~ "Upper",
        str_detect(toupper(as.character(cc50_censor)), "INTERP") ~ "Interpolated",
        is.finite(cc50) & abs(cc50 - LOW_CAP) < 1e-12 ~ "Lower",
        is.finite(cc50) & abs(cc50 - HIGH_CAP) < 1e-12 ~ "Upper",
        is.finite(cc50) ~ "Interpolated",
        TRUE ~ "QC_Failed"
      ),
      Censor_Category = factor(
        Censor_Category,
        levels=c("Lower","Interpolated","Upper","QC_Failed")
      ),
      Sensitivity_Category = factor(
        case_when(
          Censor_Category == "Lower" ~ "<=0.2",
          Censor_Category == "Interpolated" ~ "0.2-20",
          Censor_Category == "Upper" ~ ">20",
          TRUE ~ NA_character_
        ),
        levels=c("<=0.2","0.2-20",">20"),
        ordered=TRUE
      )
    )
}

primary <- read_csv(primary_file, show_col_types=FALSE) %>%
  standardize_cc50() %>%
  filter(is.finite(cc50))

strict <- read_csv(strict_file, show_col_types=FALSE) %>%
  standardize_cc50() %>%
  filter(is.finite(cc50), !str_detect(as.character(cc50_censor), "QC_FAILED"))

plate_qc <- read_csv(plate_qc_file, show_col_types=FALSE)

if (nrow(primary) != 400) warning("Expected 400 primary observations; found ", nrow(primary))
if (nrow(strict) != 272) warning("Expected 272 strict observations; found ", nrow(strict))

# ------------------------------------------------------------------------------
# 2. Descriptive summaries
# ------------------------------------------------------------------------------

summarize_cc50 <- function(dat, grouping_vars=character(0)) {
  dat %>%
    filter(is.finite(cc50)) %>%
    group_by(across(all_of(grouping_vars))) %>%
    summarise(
      n=n(),
      n_bioreps=n_distinct(BioRep_ID),
      mean_boundary=mean(cc50),
      sd_boundary=sd(cc50),
      median_boundary=median(cc50),
      q25_boundary=quantile(cc50,.25),
      q75_boundary=quantile(cc50,.75),
      n_lower=sum(Censor_Category=="Lower"),
      n_interpolated=sum(Censor_Category=="Interpolated"),
      n_upper=sum(Censor_Category=="Upper"),
      prop_lower=mean(Censor_Category=="Lower"),
      prop_interpolated=mean(Censor_Category=="Interpolated"),
      prop_upper=mean(Censor_Category=="Upper"),
      .groups="drop"
    )
}

cell_desc <- summarize_cc50(primary, "cell_line")
condition_desc <- summarize_cc50(primary, c("cell_line","Treatment","Size","Timepoint"))

write_out(cell_desc, "V6_3_CC50_ByCellLine_Descriptive.csv")
write_out(condition_desc, "V6_3_CC50_ByCondition_Descriptive.csv")

# ------------------------------------------------------------------------------
# 3. LMM helpers
# ------------------------------------------------------------------------------

fit_lmm <- function(dat, formula, label) {
  environment(formula) <- environment()
  tryCatch(
    lmerTest::lmer(formula, data=dat, REML=FALSE),
    error=function(e) {
      warning(label, " failed: ", conditionMessage(e))
      NULL
    }
  )
}

extract_type3 <- function(fit, analysis) {
  if (is.null(fit)) return(tibble())
  a <- as.data.frame(anova(fit, type=3, ddf="Satterthwaite")) %>%
    rownames_to_column("Term") %>% as_tibble()
  pcol <- grep("^Pr\\(", names(a), value=TRUE)
  fcol <- grep("^F value$|^F.value$|^F$", names(a), value=TRUE)

  a %>%
    transmute(
      Analysis=analysis,
      Term,
      NumDF=if ("NumDF" %in% names(a)) NumDF else NA_real_,
      DenDF=if ("DenDF" %in% names(a)) DenDF else NA_real_,
      F_value=if (length(fcol)) .data[[fcol[1]]] else NA_real_,
      p=if (length(pcol)) .data[[pcol[1]]] else NA_real_,
      q_BH=p.adjust(p, "BH")
    )
}

diag_lmm <- function(fit, analysis) {
  if (is.null(fit)) return(tibble())
  tibble(
    Analysis=analysis,
    N=nobs(fit),
    N_BioRep_ID=n_distinct(model.frame(fit)$BioRep_ID),
    Sigma=sigma(fit),
    Singular=lme4::isSingular(fit,tol=1e-4)
  )
}

ref_contrasts <- function(fit, dat, reference=REF_CELL, analysis) {
  if (is.null(fit)) return(tibble())

  emm <- emmeans(fit, ~ cell_line, data=dat)
  grid <- as.data.frame(emm)
  ref_idx <- match(reference, as.character(grid$cell_line))
  if (is.na(ref_idx)) stop("Reference cell line not found in emmeans grid: ", reference)

  con <- contrast(emm, method="trt.vs.ctrl", ref=ref_idx, adjust="none")

  as_tibble(summary(con, infer=c(TRUE,TRUE))) %>%
    mutate(
      Analysis=analysis,
      test_cell=sub(paste0(" - ",reference,"$"),"",contrast),
      ref_cell=reference,
      p_raw=p.value,
      q_BH=p.adjust(p_raw,"BH")
    ) %>%
    rename(
      Estimate=estimate, SE=SE, df=df, t=t.ratio,
      lowerCL=lower.CL, upperCL=upper.CL
    )
}

all_pairwise_celllines <- function(fit, dat, analysis) {
  if (is.null(fit)) return(tibble())
  emm <- emmeans(fit, ~ cell_line, data=dat)
  con <- pairs(emm, adjust="none")
  as_tibble(summary(con, infer=c(TRUE,TRUE))) %>%
    mutate(
      Analysis=analysis,
      p_raw=p.value,
      q_BH=p.adjust(p_raw,"BH")
    ) %>%
    rename(
      Estimate=estimate, SE=SE, df=df, t=t.ratio,
      lowerCL=lower.CL, upperCL=upper.CL
    )
}

factor_pairwise <- function(fit, dat, factor_name, analysis) {
  emm <- emmeans(fit, as.formula(paste0("~",factor_name)), data=dat)
  as_tibble(summary(pairs(emm,adjust="none"),infer=c(TRUE,TRUE))) %>%
    mutate(
      Analysis=analysis,
      Factor=factor_name,
      p_raw=p.value,
      q_BH=p.adjust(p_raw,"BH")
    ) %>%
    rename(
      Estimate=estimate, SE=SE, df=df, t=t.ratio,
      lowerCL=lower.CL, upperCL=upper.CL
    )
}

interaction_ref_contrasts <- function(fit, dat, by_factor, reference=REF_CELL, analysis) {
  emm <- emmeans(
    fit,
    as.formula(paste0("~ cell_line | ",by_factor)),
    data=dat
  )
  grid <- as.data.frame(emm)
  # reference index is constant because cell_line ordering repeats within panels
  ref_idx <- match(reference, levels(dat$cell_line))
  con <- contrast(emm,method="trt.vs.ctrl",ref=ref_idx,adjust="none")

  as_tibble(summary(con,infer=c(TRUE,TRUE))) %>%
    mutate(
      Analysis=analysis,
      test_cell=sub(paste0(" - ",reference,"$"),"",contrast),
      ref_cell=reference,
      p_raw=p.value
    ) %>%
    group_by(.data[[by_factor]]) %>%
    mutate(q_BH_within_panel=p.adjust(p_raw,"BH")) %>%
    ungroup() %>%
    rename(
      Estimate=estimate, SE=SE, df=df, t=t.ratio,
      lowerCL=lower.CL, upperCL=upper.CL
    )
}

# ------------------------------------------------------------------------------
# 4. Primary additive LMM
# ------------------------------------------------------------------------------

form_primary <- cc50 ~ cell_line + Treatment + Size + Timepoint + (1|BioRep_ID)
fit_primary <- fit_lmm(primary, form_primary, "Primary additive LMM")

primary_type3 <- extract_type3(fit_primary,"Primary_additive")
primary_diag <- diag_lmm(fit_primary,"Primary_additive")
primary_bewo <- ref_contrasts(fit_primary,primary,"BeWo","Primary_additive")
primary_allpairs <- all_pairwise_celllines(fit_primary,primary,"Primary_additive")
primary_treatment <- factor_pairwise(fit_primary,primary,"Treatment","Primary_additive")
primary_time <- factor_pairwise(fit_primary,primary,"Timepoint","Primary_additive")
primary_size <- factor_pairwise(fit_primary,primary,"Size","Primary_additive")

write_out(primary_type3,"V6_3_PRIMARY_Additive_LMM_TypeIII.csv")
write_out(primary_diag,"V6_3_PRIMARY_Additive_LMM_Diagnostics.csv")
write_out(primary_bewo,"V6_3_PRIMARY_CellLine_vs_BeWo.csv")
write_out(primary_allpairs,"V6_3_SUPP_AllPairwise_CellLine_Contrasts.csv")
write_out(primary_treatment,"V6_3_PRIMARY_Treatment_Pairwise.csv")
write_out(primary_time,"V6_3_PRIMARY_Timepoint_Pairwise.csv")
write_out(primary_size,"V6_3_PRIMARY_Size_Pairwise.csv")

saveRDS(fit_primary,file.path(out_dir,"V6_3_PRIMARY_Additive_LMM.rds"))

# ------------------------------------------------------------------------------
# 5. Secondary one-interaction models
# ------------------------------------------------------------------------------

fit_tx <- fit_lmm(
  primary,
  cc50 ~ cell_line*Treatment + Size + Timepoint + (1|BioRep_ID),
  "Cell line x treatment"
)
fit_time <- fit_lmm(
  primary,
  cc50 ~ cell_line*Timepoint + Size + Treatment + (1|BioRep_ID),
  "Cell line x timepoint"
)
fit_size <- fit_lmm(
  primary,
  cc50 ~ cell_line*Size + Timepoint + Treatment + (1|BioRep_ID),
  "Cell line x size"
)

tx_type3 <- extract_type3(fit_tx,"Secondary_cellline_x_Treatment")
time_type3 <- extract_type3(fit_time,"Secondary_cellline_x_Timepoint")
size_type3 <- extract_type3(fit_size,"Secondary_cellline_x_Size")

tx_bewo <- interaction_ref_contrasts(
  fit_tx,primary,"Treatment","BeWo","Secondary_cellline_x_Treatment"
)
time_bewo <- interaction_ref_contrasts(
  fit_time,primary,"Timepoint","BeWo","Secondary_cellline_x_Timepoint"
)
size_bewo <- interaction_ref_contrasts(
  fit_size,primary,"Size","BeWo","Secondary_cellline_x_Size"
)

write_out(tx_type3,"V6_3_SECONDARY_CellLine_x_Treatment_TypeIII.csv")
write_out(time_type3,"V6_3_SECONDARY_CellLine_x_Timepoint_TypeIII.csv")
write_out(size_type3,"V6_3_SECONDARY_CellLine_x_Size_TypeIII.csv")
write_out(tx_bewo,"V6_3_SECONDARY_CellLine_x_Treatment_BeWoContrasts.csv")
write_out(time_bewo,"V6_3_SECONDARY_CellLine_x_Timepoint_BeWoContrasts.csv")
write_out(size_bewo,"V6_3_SECONDARY_CellLine_x_Size_BeWoContrasts.csv")

# ------------------------------------------------------------------------------
# 6. Strict SSMD >=1 additive LMM sensitivity
# ------------------------------------------------------------------------------

fit_strict <- fit_lmm(strict,form_primary,"Strict SSMD additive LMM")
strict_type3 <- extract_type3(fit_strict,"Strict_SSMD_ge1_additive")
strict_bewo <- ref_contrasts(
  fit_strict,strict,"BeWo","Strict_SSMD_ge1_additive"
)
strict_allpairs <- all_pairwise_celllines(
  fit_strict,strict,"Strict_SSMD_ge1_additive"
)

write_out(strict_type3,"V6_3_SENS_StrictSSMD_Additive_TypeIII.csv")
write_out(strict_bewo,"V6_3_SENS_StrictSSMD_CellLine_vs_BeWo.csv")
write_out(strict_allpairs,"V6_3_SENS_StrictSSMD_AllPairwise_CellLine.csv")

# ------------------------------------------------------------------------------
# 7. Additive ordinal CLMM
# ------------------------------------------------------------------------------

ord_dat <- primary %>%
  filter(!is.na(Sensitivity_Category)) %>%
  droplevels()

ord_formula <- Sensitivity_Category ~
  cell_line + Treatment + Size + Timepoint + (1|BioRep_ID)

fit_ord <- ordinal::clmm(
  ord_formula,
  data=ord_dat,
  link="logit",
  Hess=TRUE,
  nAGQ=5
)

saveRDS(fit_ord,file.path(out_dir,"V6_3_SENS_Additive_OrdinalCLMM.rds"))

# Likelihood-ratio tests
ord_drop <- as.data.frame(drop1(fit_ord,test="Chisq")) %>%
  rownames_to_column("Term") %>%
  as_tibble()

pcol <- grep("^Pr\\(",names(ord_drop),value=TRUE)
lrtcol <- grep("^LRT$|^Chisq$|^Chisq value$",names(ord_drop),value=TRUE)

ord_lrt <- ord_drop %>%
  filter(Term!="<none>") %>%
  transmute(
    Term,
    LRT=if(length(lrtcol)) .data[[lrtcol[1]]] else NA_real_,
    p=if(length(pcol)) .data[[pcol[1]]] else NA_real_,
    q_BH=p.adjust(p,"BH")
  )

write_out(ord_lrt,"V6_3_SENS_OrdinalCLMM_LRT_Tests.csv")

# Robust ordinal emmeans helper
tidy_ord_contrast <- function(con) {
  s <- as_tibble(summary(con,infer=c(TRUE,TRUE)))
  lower_name <- intersect(c("lower.CL","asymp.LCL","lower.HPD","LCL"),names(s))
  upper_name <- intersect(c("upper.CL","asymp.UCL","upper.HPD","UCL"),names(s))

  s$CI_lower <- if(length(lower_name)) s[[lower_name[1]]] else NA_real_
  s$CI_upper <- if(length(upper_name)) s[[upper_name[1]]] else NA_real_

  s %>%
    mutate(
      p_raw=p.value,
      q_BH=p.adjust(p_raw,"BH"),
      OR_higher_CC50_category=exp(estimate),
      OR_lower95=exp(CI_lower),
      OR_upper95=exp(CI_upper)
    )
}

ord_emm_cell <- emmeans(fit_ord,~cell_line,mode="latent",data=ord_dat)
ord_grid <- as.data.frame(ord_emm_cell)
ord_ref_idx <- match("BeWo",as.character(ord_grid$cell_line))

ord_bewo <- contrast(
  ord_emm_cell,
  method="trt.vs.ctrl",
  ref=ord_ref_idx,
  adjust="none"
) %>%
  tidy_ord_contrast() %>%
  mutate(
    test_cell=sub(" - BeWo$","",contrast),
    ref_cell="BeWo"
  )

ord_allpairs <- pairs(ord_emm_cell,adjust="none") %>%
  tidy_ord_contrast()

write_out(ord_bewo,"V6_3_SENS_OrdinalCLMM_CellLine_vs_BeWo.csv")
write_out(ord_allpairs,"V6_3_SENS_OrdinalCLMM_AllPairwise_CellLine.csv")

# Polymer/time/size ordinal contrasts
for (fac in c("Treatment","Timepoint","Size")) {
  emm <- emmeans(
    fit_ord,
    as.formula(paste0("~",fac)),
    mode="latent",
    data=ord_dat
  )
  out <- pairs(emm,adjust="none") %>% tidy_ord_contrast()
  write_out(out,paste0("V6_3_SENS_OrdinalCLMM_",fac,"_Pairwise.csv"))
}

# ------------------------------------------------------------------------------
# 8. Figures
# ------------------------------------------------------------------------------

# Censoring distribution
censor_cell <- primary %>%
  count(cell_line,Censor_Category,name="n") %>%
  group_by(cell_line) %>%
  mutate(total=sum(n),proportion=n/total) %>%
  ungroup()

p_censor <- censor_cell %>%
  filter(Censor_Category!="QC_Failed") %>%
  ggplot(aes(cell_line,proportion,fill=Censor_Category)) +
  geom_col(width=.7) +
  scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
  labs(
    x=NULL,
    y="Proportion of replicate-level CC50 estimates",
    fill="CC50 category",
    title="Assay-range censoring by cell line"
  ) +
  theme_classic(base_size=12)

ggsave(
  file.path(plot_dir,"Fig_V6_3_CC50_Censoring_ByCellLine.png"),
  p_censor,width=7,height=4.8,dpi=400,bg="white"
)

# Distribution plot
p_dist <- ggplot(primary,aes(cell_line,cc50)) +
  geom_boxplot(outlier.shape=NA,width=.6) +
  geom_jitter(width=.13,height=0,alpha=.4,size=1.2) +
  geom_hline(yintercept=c(LOW_CAP,HIGH_CAP),linetype="dotted") +
  scale_y_continuous(breaks=c(.2,2,5,10,15,20)) +
  coord_cartesian(ylim=c(0,HIGH_CAP)) +
  labs(
    x=NULL,
    y=expression(CC[50]~"(µg/mL; boundary-assigned)"),
    title="QC-valid replicate-level CC50 estimates by cell line",
    subtitle="Values <=0.2 and >20 µg/mL are displayed at assay boundaries"
  ) +
  theme_classic(base_size=12)

ggsave(
  file.path(plot_dir,"Fig_V6_3_CC50_Distribution_ByCellLine.png"),
  p_dist,width=7.5,height=5,dpi=400,bg="white"
)

# BeWo contrast plot
p_bewo <- primary_bewo %>%
  mutate(sig=q_BH<ALPHA,label=paste0(test_cell," vs ",ref_cell)) %>%
  ggplot(aes(x=Estimate,y=label)) +
  geom_vline(xintercept=0,linetype="dashed") +
  geom_errorbar(aes(xmin=lowerCL,xmax=upperCL),orientation="y",width=.2) +
  geom_point(aes(shape=sig),size=2.6) +
  labs(
    x=expression(Delta*CC[50]~"vs BeWo (µg/mL)"),
    y=NULL,
    shape="BH q < 0.05",
    title="Adjusted cell-line differences versus BeWo",
    subtitle="Primary additive mixed-effects model"
  ) +
  theme_classic(base_size=12)

ggsave(
  file.path(plot_dir,"Fig_V6_3_PRIMARY_CellLine_vs_BeWo.png"),
  p_bewo,width=7.2,height=4.5,dpi=400,bg="white"
)

# ------------------------------------------------------------------------------
# 9. Results-ready summary
# ------------------------------------------------------------------------------

fmtp <- function(x) {
  ifelse(is.na(x),"NA",ifelse(x<.001,"<0.001",sprintf("%.3f",x)))
}

lines <- c(
  "MNP CC50 v6.3 FINAL — BeWo reference",
  "===================================",
  paste0("Primary observations: ",nrow(primary)),
  paste0("Strict SSMD>=1 observations: ",nrow(strict)),
  "",
  "PRIMARY MODEL:",
  "cc50 ~ cell_line + Treatment + Size + Timepoint + (1|BioRep_ID)",
  "",
  "Primary Type III tests:"
)

for(i in seq_len(nrow(primary_type3))) {
  r <- primary_type3[i,]
  lines <- c(lines,sprintf(
    "  %s: F=%s, p=%s, BH q=%s",
    r$Term,
    ifelse(is.na(r$F_value),"NA",sprintf("%.3f",r$F_value)),
    fmtp(r$p),fmtp(r$q_BH)
  ))
}

lines <- c(lines,"","Primary adjusted contrasts versus BeWo:")
for(i in seq_len(nrow(primary_bewo))) {
  r <- primary_bewo[i,]
  lines <- c(lines,sprintf(
    "  %s vs BeWo: estimate %.2f (95%% CI %.2f to %.2f), p=%s, BH q=%s",
    r$test_cell,r$Estimate,r$lowerCL,r$upperCL,fmtp(r$p_raw),fmtp(r$q_BH)
  ))
}

lines <- c(lines,"","Ordinal sensitivity LRT tests:")
for(i in seq_len(nrow(ord_lrt))) {
  r <- ord_lrt[i,]
  lines <- c(lines,sprintf(
    "  %s: LRT=%s, p=%s, BH q=%s",
    r$Term,
    ifelse(is.na(r$LRT),"NA",sprintf("%.3f",r$LRT)),
    fmtp(r$p),fmtp(r$q_BH)
  ))
}

lines <- c(lines,"","Ordinal adjusted contrasts versus BeWo:")
for(i in seq_len(nrow(ord_bewo))) {
  r <- ord_bewo[i,]
  lines <- c(lines,sprintf(
    "  %s vs BeWo: OR %.3f (95%% CI %.3f to %.3f), p=%s, BH q=%s",
    r$test_cell,
    r$OR_higher_CC50_category,
    r$OR_lower95,
    r$OR_upper95,
    fmtp(r$p_raw),
    fmtp(r$q_BH)
  ))
}

writeLines(lines,file.path(out_dir,"V6_3_RESULTS_READY_SUMMARY.txt"))

capture.output(sessionInfo(),file=file.path(out_dir,"V6_3_sessionInfo.txt"))

cat("\n==================== v6.3 COMPLETE ====================\n")
cat("Output directory: ",out_dir,"\n",sep="")
cat("Reference cell line: BeWo\n")
cat("Primary observations: ",nrow(primary),"\n",sep="")
cat("Strict SSMD>=1 observations: ",nrow(strict),"\n",sep="")
cat("Primary LMM singular: ",lme4::isSingular(fit_primary,tol=1e-4),"\n",sep="")
cat("Treatment-interaction singular: ",lme4::isSingular(fit_tx,tol=1e-4),"\n",sep="")
cat("Timepoint-interaction singular: ",lme4::isSingular(fit_time,tol=1e-4),"\n",sep="")
cat("Size-interaction singular: ",lme4::isSingular(fit_size,tol=1e-4),"\n",sep="")
cat("Ordinal additive model fitted: TRUE\n")
cat("All pairwise cell-line contrasts exported: TRUE\n")
cat("=======================================================\n")
