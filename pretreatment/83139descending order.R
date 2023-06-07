rm(list = ls())
install.packages(dplyr)
install.packages(Seurat)
install.packages(patchwork)
#install.packages("readxl")
library(dplyr)
library(Seurat)
library(patchwork)
#library(readxl)
data<-read.csv('E:\\毕业论文\\单细胞数据集\\GSE83139\\GSE83139.csv',,header=T,row.names=1,encoding = "UTF-8")#第一行为列名,第一列为行名
sxjybl<-CreateSeuratObject(counts =data,project = "data", min.cells = 635, min.features = 200, names.delim = "_",)
ncol(sxjybl)#
nrow(sxjybl)#
A<-sxjybl@assays$RNA@counts[1:19950,1:635]
b<-as.data.frame(A)
expres=b[order(rowSums(b),decreasing = T),]
write.csv(expres,file ="E:\\毕业论文\\单细胞数据集\\GSE83139\\GSE83139（降序后）.csv")
