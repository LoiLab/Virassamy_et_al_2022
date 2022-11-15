library(Seurat)
library(reshape2)
library(ruvIIInb)
library(SingleCellExperiment)
library(scater)
library(scran)
library(tidyverse)
library(bluster)
library(readxl)
library(GSEABase)
library(AUCell)
library(ggrepel)

load("~/tils.obj.rdata")


# AUCell score
tils.count <- GetAssayData(tils.obj, assay = "RNA", slot = "data")
tils.count <- as.matrix(tils.count)
signatures <- read_excel("~/Supplementary1.xlsx")
trm_w6 <- signatures$Symbol...1
trm_w6 <- trm_w6[!is.na(trm_w6)]
tex_w6 <- signatures$Symbol...4
w6.set <- list()
w6.set[[1]] <- GeneSet(trm_w6, setName = "TRM")
w6.set[[2]] <- GeneSet(tex_w6, setName = "TEX")
w6.set <- GeneSetCollection(w6.set)
w6.auc <- AUCell_run(tils.count, w6.set)

# Wilcoxon rank sum test
trm.itgae.w6s <- getAUC(w6.auc)["TRM",which(tils.obj$anno %in% "CD8+ITGAE+")]
trm.tox.w6s <- getAUC(w6.auc)["TRM",which(tils.obj$anno %in% "CD8+TOXhigh")]
trm.w6p <- wilcox.test(trm.itgae.w6s, trm.tox.w6s)
tex.itgae.w6s <- getAUC(w6.auc)["TEX",which(tils.obj$anno %in% "CD8+ITGAE+")]
tex.tox.w6s <- getAUC(w6.auc)["TEX",which(tils.obj$anno %in% "CD8+TOXhigh")]
tex.w6p <- wilcox.test(tex.itgae.w6s, tex.tox.w6s)


# Jaccard index
vehicle.clone.anno <- read.csv("~/all_contig_annotations.csv")
cell_clone_info <- function(obj, clone.anno){
  clone.info <- clone.anno[,c(1,6,7,9,13,15,16,17,18)]
  clone.info <- clone.info[which(clone.info$chain %in% "TRB"),]
  # Filter clones with none consensus id
  clone.info <- clone.info[which(!clone.info$raw_consensus_id %in% "None"),]
  # Select each clone with largest reads
  clone.info <- clone.info %>% 
    group_by(barcode, raw_clonotype_id) %>% 
    top_n(1, reads)
  clone.info <- clone.info[,c(1,3,4,5,8)]
  clone.info$barcode <- gsub("-1", "", clone.info$barcode)
  colnames(clone.info)[5] <- "clonotype" 
  clone.info <- clone.info[which(clone.info$barcode %in% colnames(obj)),]
  clusters <- obj$anno
  names(clusters) <- colnames(obj)
  clone.info$cluster <- clusters[clone.info$barcode]
  return(clone.info)
}
tils.clone1 <- cell_clone_info(tils.obj[,which(tils.obj$sample %in% "mouse1")], vehicle.clone.anno)
tils.clone2 <- cell_clone_info(tils.obj[,which(tils.obj$sample %in% "mouse2")], vehicle.clone.anno)
tils.clone3 <- cell_clone_info(tils.obj[,which(tils.obj$sample %in% "mouse3")], vehicle.clone.anno)

VDJOverlap <- function(clone){
  cdr.freq.cl <- clone %>% 
    group_by(cdr3, cluster) %>% 
    summarise(freq = n(), .groups = 'drop')
  all.com <- cdr.freq.cl %>% 
    tidyr::expand(cdr3, cluster) 
  cdr.freq.cl2 <- merge(cdr.freq.cl, all.com, by = c("cdr3", "cluster"), all.x = TRUE, all.y = TRUE)
  cdr.freq.cl2$freq[is.na(cdr.freq.cl2$freq)] <- 0
  data.mat <- matrix(NA, nrow = length(levels(cdr.freq.cl2$cluster)), ncol = length(levels(cdr.freq.cl2$cluster)))
  data.list <- list("MH" = data.mat, "Jaccard" = data.mat)
  for(i in 1:(nrow(data.mat)-1)){ 
    for(j in (i+1):ncol(data.mat)){ 
      sub.x.df <- cdr.freq.cl2[which(cdr.freq.cl2$cluster %in% levels(cdr.freq.cl2$cluster)[i]),]
      sub.y.df <- cdr.freq.cl2[which(cdr.freq.cl2$cluster %in% levels(cdr.freq.cl2$cluster)[j]),]
      sub.x <- sub.x.df$freq
      sub.y <- sub.y.df$freq
      X <- sum(sub.x)
      Y <- sum(sub.y)
      sub.x.df2 <- sub.x.df[which(sub.x.df$freq > 0),]
      sub.y.df2 <- sub.y.df[which(sub.y.df$freq > 0),]
      cdr3_inter <- intersect(sub.x.df2$cdr3, sub.y.df2$cdr3)
      if(length(cdr3_inter) > 0){ 
        sub.x.df2 <- sub.x.df2[match(cdr3_inter, sub.x.df2$cdr3),]
        sub.y.df2 <- sub.y.df2[match(cdr3_inter, sub.y.df2$cdr3),]
        x_y_inter <- sum(pmin(sub.x.df2$freq, sub.y.df2$freq))
        x_y_union <- sum(sub.x.df2$freq, sub.y.df2$freq)
        jd.idx <- x_y_inter/(x_y_union-x_y_inter)
      }else{ 
        jd.idx <- 0
      }
      data.list[["Jaccard"]][i,j] <- jd.idx
    }
  }
  for(i in 1:length(data.list)){ 
    rownames(data.list[[i]]) <- levels(cdr.freq.cl2$cluster)
    colnames(data.list[[i]]) <- levels(cdr.freq.cl2$cluster)
    data.list[[i]] <- data.list[[i]][1:(nrow(data.list[[i]])-1),2:ncol(data.list[[i]])]
  }
  return(data.list)
}
clone.overlap1 <- VDJOverlap(tils.clone1)
clone.overlap2 <- VDJOverlap(tils.clone2)
clone.overlap3 <- VDJOverlap(tils.clone3)







