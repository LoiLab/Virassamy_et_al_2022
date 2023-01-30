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


# Liu et al ---------------------------------------------------------------

# The data required for these calculations is available in the supplementary tables of "Liu D, Schilling B, Liu D, Sucker A, Livingstone E, Jerby-Amon L, Zimmer L, Gutzmer R, Satzger I, Loquai C, Grabbe S, Vokes N, Margolis CA, Conway J, He MX, Elmarakeby H, Dietlein F, Miao D, Tracy A, Gogas H, Goldinger SM, Utikal J, Blank CU, Rauschenberg R, von Bubnoff D, Krackhardt A, Weide B, Haferkamp S, Kiecker F, Izar B, Garraway L, Regev A, Flaherty K, Paschen A, Van Allen EM, Schadendorf D (2019) Integrative molecular and clinical modeling of clinical outcomes to PD1 blockade in patients with metastatic melanoma. Nat Med 25:1916–1927. doi: 10.1038/s41591-019-0654-5"
# The clinical data is the first worksheet in the 'Supplementary Tables' Excel file (41591_2019_654_MOESM4_ESM.xlsx)

clinical_data.df <- read_excel(path = "41591_2019_654_MOESM4_ESM.xlsx", sheet = 1, skip = 2) %>%
  filter(!is.na(Tx)) %>%
  dplyr::rename(case = `...1`)

# The RNASeq gene counts (TPM format) are in the 'Supplementary Data 2' tab separated value text file (41591_2019_654_MOESM3_ESM.txt)

gexp.df <- read_tsv(file = "41591_2019_654_MOESM3_ESM.txt") %>%
  dplyr::rename(case = `...1`)

# prepare dataset
liu_clin_sig.df <- gexp.df %>%
  dplyr::select(case, one_of(sig.df$Symbol)) %>%
  gather(key = gene, value = tpm, -case) %>%
  group_by(case) %>%
  # generate signature score per case, calculated as mean of signature gene TPM. Equivalent to genefu 'sig.score' function
  dplyr::summarise(mean_exp = mean(tpm)) %>%
  left_join(clinical_data.df,.) %>%
  # select cases with prior anti-CTLA4 therapy
  filter(priorCTLA4 == 1) %>%
  filter(!is.na(mean_exp)) %>%
  # recode best response for plotting
  mutate(BR2 = ifelse(BR %in% c("PR","CR"), 1, 0)) %>%
  dplyr::select(case, PFS, progressed, OS, dead, BR2, mean_exp) %>%
  # dichotomise by median signature score
  mutate(sig_cat = ifelse(sig < median(sig), "TRMSig Low","TRMSig High"))



#### ROC curve ####
proc <- roc(liu_clin_sig.df$BR2, liu_clin_sig.df$mean_exp,
            smoothed = T,
            ci = T, ci.alpha = 0.9, stratified = F,
            plot = T, auc.polygon = T, max.auc.polygon = T, grid = T, print.auc = T, show.thres = T)

proc_auc <- round(as.numeric(proc$ci),2)

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
       title = "D40 TRMsig and ORR with anti-PD1 post\nanti-CTLA4 in metastatic melanoma") +
  annotate(geom = "text", x = 1, y = 0.05, label = paste("n = ", proc_ncases, "\nAUC ", proc_auc[2], collapse = "", sep = ""), size = 5, hjust = 1) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 16),
        panel.grid = element_blank(),
        plot.title = element_text(size = 16, hjust = 0.5))



#### Survival Analysis ####
km_liu_pfs <- survfit(Surv(PFS/(365.24/12),progressed) ~ sig_cat, data = liu_clin_sig.df) %>%
  ggsurvplot(pval = F,
             risk.table = T,
             censor.shape = 124,
             palette = c("#8B1A4F","#727C34"),
             size = 1.5,
             legend = "none",
             title = "D40 TRMsig and PFS with anti-PD1 post anti-CTLA4\nin metastatic melanoma",
             xlab = "Time (months)",
             ylab = "Progression-free Survival",
             break.time.by = 12)


# extract summary statistics (logrank test)
km_liu_pfs.tidy <- broom::tidy(survdiff(Surv(PFS/(365.24/12),progressed) ~ sig_cat, data = liu_clin_sig.df))
km_liu_pfs.glance <- broom::glance(survdiff(Surv(PFS/(365.24/12),progressed) ~ sig_cat, data = liu_clin_sig.df))

# produce final plot for publication
km_liu_pfs$plot <- ggpar(km_liu_pfs$plot, font.title = 18, font.tickslab = 16, font.x = 16, font.y = 16) + 
  annotate(geom = "text", x = 0, y = 0.05, label = sprintf("n = %s, p = %s", sum(km_liu_pfs.tidy$N),signif(km_liu_pfs.glance$p.value, digits = 2)), hjust = 0, size = 5) +
  annotate(geom = "text", x = 32, y = 1.0, label = "D40 TRM Signature high", hjust = 0, size = 5) +
  annotate(geom = "segment", x = 29, y = 1.0, xend = 31, yend = 1.0, size = 1.5, color = "#8B1A4F" ) +
  annotate(geom = "text", x = 32, y = 0.9, label = "D40 TRM Signature low", hjust = 0, size = 5) +
  annotate(geom = "segment", x = 29, y = 0.9, xend = 31, yend = 0.9, size = 1.5, color = "#727C34" ) +
  theme(plot.title = element_text(hjust = 0.5))

km_liu_pfs$plot
