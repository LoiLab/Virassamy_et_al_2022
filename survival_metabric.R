rm(list = ls())
library(dplyr)
library(data.table)
library(genefu)
library(ggplot2)

## Human homologs from Table S5, some symbols have changed in the metabric microarray annotation
signature_genes = c("FLJ14166","FOXP3","PID1","SIGLEC1","MRC1","ITGAE","GZMA","CYP4F3","RETNLB","CD300LF","F13A1","LPL","PLTP","RASSF4","TSPAN9","FCGR2B","FCGR2A","FCGR2C","SLC11A1","CMKLR1","CCL23","CCL15","SIRPA","SIRPG","SIRPB1","DAB2","LY86","GZMB","CSF1R","ITGAM","STAB1","TGFBI","FLJ22662","PLD4","FES","ST8SIA1","MS4A7","MAFB","C3AR1","CENTG2","LGMN","SASH1","C1QB","CD74","SYK","APOE","HLA-DQB1","HLA-DQB2","ITSN1","APOL2","APOL1","APOL4","APOL3","CPD","SRC","ITGB7")

## Combined expression of signature genes and clinical data from CBioportal
metabric_tbl = fread('~/Documents/Work/repositories/Virassamy_et_al_2022/metabric.csv') %>% data.frame(check.names = F) 

## Calculate signature scores
x = data.frame(probe = signature_genes,EntrezGene.ID=signature_genes,coefficient=1)
annot = data.frame(EntrezGene.ID=signature_genes)
scores = sig.score(x = x, data = metabric_tbl[,signature_genes], annot = annot)
metabric_tbl$signature_score = scores$score

## Survival Plots

# Filter Basal samples
metabric_tbl = metabric_tbl %>% filter(Pam50Subtype=='Basal')
metabric_tbl$signature_score_d = ifelse(metabric_tbl$signature_score>median(metabric_tbl$signature_score),'high','low')


## DDFS
sfit <- survfit(Surv(t.ddfs, e.ddfs) ~ signature_score_d, data=metabric_tbl)
coxregre <- coxph(Surv(t.ddfs,e.ddfs)~signature_score,data = metabric_tbl)
hr = signif(summary(coxregre)$coef[1,"exp(coef)"],digits=2)
ci.up= signif(summary(coxregre)$conf.int[1,4],digits=2)
ci.down = signif(summary(coxregre)$conf.int[1,3],digits=2)
pval = signif(summary(coxregre)$coefficients[1,5],digits=4)

hr.txt = c('HR: 95% CI')
pval.txt = c('p value')
hr.txt = paste(hr.txt,paste(hr,": ",ci.down,"-",ci.up,sep = ""))
pval.txt = paste(pval.txt,pval)

ggsurvplot(data = metabric_tbl,
                  sfit,                     # survfit object with calculated statistics.
                  conf.int = F,
                  xlab = "Time (Years)",
                  font.x = c(8, "bold", "black"),
                  ylab = "Distant Disease Free Survival",
                  linetype = c(1,2),    # customize X axis label.
                  break.time.by = 1,     # break X axis in time intervals by 200.
                  ggtheme = theme_light(), # customize plot and risk table with a theme.
                  tables.theme = theme_cleantable(),
                  fontsize = 3,
                  risk.table.y.text.col = T,# colour risk table text annotations.
                  risk.table.y.text = FALSE,# show bars instead of names in text annotations
                  title = "Metabric: DDFS Basal D40 Trm signature",
                  subtitle = paste(hr.txt,pval.txt,sep = ' '),
                  font.subtitle = c(8, "bold.italic", "purple"),
                  font.title = c(10, "bold", "blue"),
                  surv.median.line = "none",  # add the median survival pointer.
                  xlim = c(0, 10),
                 

)

## OS
sfit <- survfit(Surv(t.os, e.os) ~ signature_score_d, data=metabric_tbl)
coxregre <- coxph(Surv(t.os,e.os)~signature_score,data = metabric_tbl)
hr = signif(summary(coxregre)$coef[1,"exp(coef)"],digits=2)
ci.up= signif(summary(coxregre)$conf.int[1,4],digits=2)
ci.down = signif(summary(coxregre)$conf.int[1,3],digits=2)
pval = signif(summary(coxregre)$coefficients[1,5],digits=4)

hr.txt = c('HR: 95% CI')
pval.txt = c('p value')
hr.txt = paste(hr.txt,paste(hr,": ",ci.down,"-",ci.up,sep = ""))
pval.txt = paste(pval.txt,pval)

ggsurvplot(data = metabric_tbl,
           sfit,                     # survfit object with calculated statistics.
           conf.int = F,
           xlab = "Time (Years)",
           font.x = c(8, "bold", "black"),
           ylab = "Overall Survival",
           linetype = c(1,2),    # customize X axis label.
           break.time.by = 1,     # break X axis in time intervals by 200.
           ggtheme = theme_light(), # customize plot and risk table with a theme.
           tables.theme = theme_cleantable(),
           fontsize = 3,
           risk.table.y.text.col = T,# colour risk table text annotations.
           risk.table.y.text = FALSE,# show bars instead of names in text annotations
           title = "Metabric: OS Basal D40 Trm signature",
           subtitle = paste(hr.txt,pval.txt,sep = ' '),
           font.subtitle = c(8, "bold.italic", "purple"),
           font.title = c(10, "bold", "blue"),
           surv.median.line = "none",  # add the median survival pointer.
           xlim = c(0, 10),
           
           
)
