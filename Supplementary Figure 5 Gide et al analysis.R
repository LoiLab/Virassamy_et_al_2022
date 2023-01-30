library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(readxl)
library(pROC)
library(survival)
library(survminer)

# signature from Supplementary Table 5
sig.df <- read_excel(sheet = "Table S5", path = "Table S5.xlsx", range = "R5:S62")


# Gide et al ---------------------------------------------------------------

# The clinical data required for these calculations is available in Supplementary Table S2 (mmc3.xlsx) from: "Gide TN, Quek C, Menzies AM, Tasker AT, Shang P, Holst J, Madore J, Lim SY, Velickovic R, Wongchenko M, Yan Y, Lo S, Carlino MS, Guminski A, Saw RPM, Pang A, McGuire HM, Palendira U, Thompson JF, Rizos H, Silva IP da, Batten M, Scolyer RA, Long GV, Wilmott JS (2019) Distinct Immune Cell Populations Define Response to Anti-PD-1 Monotherapy and Anti-PD-1/Anti-CTLA-4 Combined Therapy. Cancer Cell 35:238–255.e6. doi: 10.1016/j.ccell.2019.01.003"
# this table contains data for combination Ipilimumab and anti-PD1 therapy only.
gide_clin.df <- read_excel(path = "mmc3.xlsx", skip = 2) %>%
  mutate(sample_name = ifelse(grepl("PRE", `RNA Sequencing`, fixed = T) | `Patient no.` == 38, paste("ipiPD1", `Patient no.`, "PRE", sep = "_"), "-"))


# The authors do not provide processed RNASeq data. The raw data is available at https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB23709
# The methods applied to processing this data are described in the paper methods section
# The output of the above is a matrix of gene counts, with genes are row names and sample names in columns. This is referred to as 'counts.mat'
# Note the sample names in Bio Project of the form 'PD1_NN_PRE' where NN maps to the 'Patient no.' column in the clinical data file. The suffix PRE or EDT denotes whether the sample was taken before or at the end of treatment. Only baseline samples in the combination therapy cohort were used, and these all have sample names of the form "ipiPD1_NN_PRE'.


# chose signature genes available in count matrix
sig_genes <- sig.df$Symbol[sig.df$Symbol %in% rownames(counts.mat)]

# prepare dataset
gide_clin_sig.df <- counts.mat[sig_genes,] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  gather(key = sample_name, value = exp, -gene) %>%
  group_by(gene) %>%
  mutate(exp = scale(exp, center = T, scale = T)) %>%
  ungroup %>%
  group_by(sample_name) %>%
  # generate signature score per case, calculated as mean of signature gene TPM. Equivalent to genefu 'sig.score' function
  summarise(sig = mean(exp)) %>%
  left_join(gide_clin.df) %>%
  mutate(Response = ifelse(`Best RECIST response` %in% c("SD","PD"), 0, 1)) %>%
  # dichotomise by median signature score
  mutate(sig_cat = ifelse(sig < median(sig), "TRMSig Low","TRMSig High")) %>%
  mutate(PD = ifelse(Progressed == "Yes", 1, 0)) %>%
  dplyr::rename(PFS = `Progression Free Survival (Days)`)


#### ROC curve ####
proc <- roc(gide_clin_sig.df $Response, gide_clin_sig.df$sig,
            smoothed = T,
            ci = T, ci.alpha = 0.9, stratified = F,
            plot = T, auc.polygon = T, max.auc.polygon = T, grid = T, print.auc = T, show.thres = T)


proc_auc <- round(as.numeric(proc$ci), 2)

proc_ncases = length(proc$cases) + length(proc$controls)

# function to extract data for plotting ROC curves
roc2df <- function(p) {
  
  tibble(TPR = p$sensitivities,
         FPR = 1 - p$specificities)
  
}


roc2df(proc) %>%
  distinct %>%
  dplyr::arrange(TPR) %>%
  ggplot(aes(x = FPR, y = TPR)) + 
  geom_line(size = 1.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 1, size = 1, colour = "grey50") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  labs(x = "1 - Specificity", y = "Sensitivity",
       title = "D40 TRMsig and ORR with anti-CTLA4 +\nanti-PD1 in metastatic melanoma") +
  annotate(geom = "text", x = 1, y = 0.05, label = paste("n = ", proc_ncases, "\nAUC ", proc_auc[2], collapse = "", sep = ""), size = 5, hjust = 1) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 16),
        panel.grid = element_blank(),
        plot.title = element_text(size = 16, hjust = 0.5))


#### Survival Analysis ####
km_gide_pfs <- ggsurvplot(survfit(Surv(PFS/(365.24/12), event = PD) ~ sig_cat, data = gide_clin_sig.df),
                          risk.table = T,
                          pval = F,
                          size = 1.5,
                          censor.shape = 124,
                          palette = c("#8B1A4F","#727C34"),
                          legend = "none",
                          title = "D40 TRMsig and PFS with anti-CTLA4 + anti-PD1 in\nmetastatic melanoma",
                          xlab = "Time (months)",
                          break.time.by = 12)

# extract summary statistics (logrank test)
km_gide_pfs.tidy <- broom::tidy(survdiff(Surv(PFS/(365.24/12), PD) ~ sig_cat, data = gide_clin_sig.df))
km_gide_pfs.glance <- broom::glance(survdiff(Surv(PFS/(365.24/12), PD) ~ sig_cat, data = gide_clin_sig.df))

# produce final plot for publication
km_gide_pfs$plot <- ggpar(km_gide_pfs$plot, font.title = 18, font.tickslab = 16, font.x = 16, font.y = 16) + 
  annotate(geom = "text", x = 0, y = 0.05, label = sprintf("n = %s, p = %s", sum(km_gide_pfs.tidy$N),signif(km_gide_pfs.glance$p.value, digits = 2)), hjust = 0, size = 5) +
  annotate(geom = "text", x = 17, y = 1.0, label = "D40 TRM Signature high", hjust = 0, size = 5) +
  annotate(geom = "segment", x = 14, y = 1.0, xend = 16, yend = 1.0, size = 1.5, color = "#8B1A4F" ) +
  annotate(geom = "text", x = 17, y = 0.9, label = "D40 TRM Signature low", hjust = 0, size = 5) +
  annotate(geom = "segment", x = 14, y = 0.9, xend = 16, yend = 0.9, size = 1.5, color = "#727C34" ) +
  theme(plot.title = element_text(hjust = 0.5))

km_gide_pfs$plot

