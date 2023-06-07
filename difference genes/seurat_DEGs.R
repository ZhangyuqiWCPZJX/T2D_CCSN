rm(list=ls())
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("PFP")
BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")
setwd("E:\\毕业论文\\单细胞数据集\\GSE83139\\差异基因\\beta\\seurat")
library(Seurat)

options(stringsAsFactors = F)

raw.data1<-read.csv("E:\\毕业论文\\单细胞数据集\\GSE83139\\网络构建\\度矩阵\\beta\\GSE83139beta-child-cndm.csv",header=T,row.names = 1)
raw.data2<-read.csv("E:\\毕业论文\\单细胞数据集\\GSE83139\\网络构建\\度矩阵\\beta\\GSE83139beta-T2D-cndm.csv",header=T,row.names = 1)

datan1 = data.frame(raw.data1)
datan2 = data.frame(raw.data2)

dataan1 <- as(as.matrix(datan1), "dgCMatrix")
dataan2 <- as(as.matrix(datan2), "dgCMatrix")

child<- CreateSeuratObject(counts = dataan1)#可筛选细胞和基因数
T2D<- CreateSeuratObject(counts = dataan2)

Idents(child) <- "child"
Idents(T2D) <- "T2D"

scRNA <- merge(child,T2D)
list_data <- list()
list_data[["mait"]] <- scRNA

library(PFP)
library(clusterProfiler)
library(org.Hs.eg.db)
list_diff_genes2 <- list()
i=1
for (data0 in list_data){
  DefaultAssay(data0) <- "RNA"  #设置为integrated是批处理调整的数据，设置为RNA是基于原始值分析
  #Idents(data0) <- "group"
  #genes1 <- FindMarkers(object = data0,ident.1 = "HD",ident.2 = "Moderate",logfc.threshold = 0.25)
  genes1 <- FindMarkers(object = data0,ident.1 = "child",ident.2 = "T2D" ,logfc.threshold = 0.25)
  genes1 <- genes1[order(genes1$avg_log2FC,decreasing = T),]
  genes1 <- genes1[genes1$p_val_adj<0.05,]
  
  genes2 <- bitr(geneID = rownames(genes1),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  list_diff_genes2[[i]] <- genes2$ENTREZID
  i=i+1
}

write.csv(genes1,file="beta(child-T2D)du.csv")
#write.csv(genes2,file="child-T2Ddutype.csv")
