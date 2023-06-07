rm(list=ls())
#BiocManager::install('Seurat')
library(Seurat)
library(dplyr)
library(patchwork)
# library(scCATCH)
BiocManager::install('profvis')
library(devtools)
BiocManager::install('SingleR')
#devtools::install_github('dviraran/SingleR')
library(SingleR)
library(mindr)
library(Matrix)
library(magrittr)
library(hdf5r)
library(SeuratObject)
# Load data 
#setwd("E:\\毕业论文\\单细胞数据集\\GSE83139\\聚类\\基因表达矩阵\\resolution=0.8(8)\\最终")
#setwd("E:\\毕业论文\\单细胞数据集\\GSE83139\\聚类\\基因表达矩阵")
setwd("E:\\毕业论文\\单细胞数据集\\GSE83139\\聚类\\度矩阵\\resolution=0.8（8）")
#raw.data<-read.csv("E:\\毕业论文\\单细胞数据集\\GSE83139\\网络构建\\基因表达矩阵\\GSE83139（1.5w）.csv",header=T,row.names = 1)
raw.data<-read.csv("E:\\毕业论文\\单细胞数据集\\GSE83139\\网络构建\\度矩阵\\GSE83139-cndm.csv",header=T,row.names = 1)
datan = data.frame(raw.data)
dataan <- as(as.matrix(datan), "dgCMatrix")
pbmc<- CreateSeuratObject(counts = dataan)#可筛选细胞和基因数
dim(pbmc)
## =============5.鉴定高变基因
# 高变基因：在一些细胞中表达高，另一些细胞中表达低的基因
# 变异指标： mean-variance relationship
# 默认返回两千个高变基因，用于下游如PCA降维分析。
pbmc<-na.omit(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# 提取前10的高变基因
#要对基因有一定的熟悉度，可能是某些细胞的marker
top10 <- head(VariableFeatures(pbmc), 10)
top10
top30 <- head(VariableFeatures(pbmc), 30)
top30
# 展示高变基因
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)#有问题
plot1 + plot2
#横轴代表均值，纵轴代表方差，红色前2000个高变基因



## =============6.Scaling the data
# 归一化处理：每一个基因在所有细胞中的均值变为0，方差标为1，对于降维来说是必需步骤
# 归一化后的值保存在：pbmc[["RNA"]]@scale.data
#ScaleData 默认只对2000个基因进行归一化
pbmc <- ScaleData(pbmc)
scale.data <- pbmc[["RNA"]]@scale.data
dim(scale.data)
scale.data[1:10,1:4]


# 可以选择全部基因归一化
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)




##7.降维
# PCA降维，默认使用前面2000个高变基因，可以使用features改变用于降维的基因集
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# 可视化
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# ## =============8.确定使用PC个数
# # each PC essentially representing a ‘metafeature’
# pbmc <- JackStraw(pbmc, num.replicate = 300)
# pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
# #每条线代表一个PC，颜色代表不同的Pvalue值
# JackStrawPlot(pbmc, dims = 1:5)
#贡献度
ElbowPlot(pbmc)

## =============9.对细胞聚类
# 首先基于PCA空间构建一个基于欧氏距离的KNN图，然后基于其局部邻域中的共享重叠来细化任意两个细胞之间边缘的权重（Jaccard相似性）
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# 聚类并最优化
# resolution参数：值越大，细胞分群数越多，
# 其中resolution的参数可以自己调节，resolution越大，所聚出的类的数目也越多，
#一般来说对于3k个细胞，推荐值为0.4-1.2，可根据自己的数据和先验知识或直觉进行调整
pbmc <- FindClusters(pbmc, resolution =0.8)

# 查看聚类数ID
head(Idents(pbmc), 5)

set.seed(123)
## =============10.将细胞在低维空间可视化UMAP/tSNE
pbmc <- RunUMAP(pbmc, dims = 1:10)#dims=1:10 使用PCA的某些维度作为UMAP的输入，还可以设置 features=someGenes 使用这些基因作为UMAP的输入
pbmc <- RunTSNE(pbmc, dims = 1:10, check_duplicates = FALSE)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)#比较一个cluster与所有其他cluster之间的基因表达（差异基因）
write.csv(pbmc.markers,file = "pbmc.markers.csv")
library(ggplot2)
library(dplyr) 
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
p0<-DoHeatmap(pbmc,top10$gene,size=5)
# DoHeatmap(pbmc,top10$gene,size=5,lines.width = 4,group.bar.height = 0.1)
ggsave("pbmc.markers各类表达热图.png",p0,width=15,height=45)

# 可视化，使用UMAP多一点

p1=DimPlot(pbmc, reduction = "umap", label = T, label.size = 5)
ggsave("CellType.UMAP.png",p1,width=10,height=8)
p2=DimPlot(pbmc, reduction = "tsne", label = T, label.size = 5)
ggsave("CellType.TSNE.png",p2,width=10,height=8)

#保存pbmc
saveRDS(pbmc, file = "pbmc_du.rds")#需要改
#pbmc<-readRDS("E:\\毕业论文\\单细胞数据集\\GSE83139\\聚类\\基因表达矩阵\\resolution=0.8(8)\\pbmc.rds")

##==鉴定细胞类（自动注释）型==##=======================================================================
library(SingleR)
# BiocManager::install('celldex')
library(celldex)
#refdata <-HumanPrimaryCellAtlasData()
#refdata <- MonacoImmuneData()
#refdata <- BlueprintEncodeData()
#refdata <- ImmGenData()
#refdata <- DatabaseImmuneCellExpressionData()
refdata <- NovershternHematopoieticData()
testdata <- GetAssayData(pbmc, slot="data")  #相当于pbmc@assays$RNA@data
# pbmc_SingleR <- pbmc@assays$RNA@counts
clusters <- pbmc@meta.data$seurat_clusters
#使用HPCA参考数据库鉴定
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"celltype_HPCA.csv",row.names = F)
saveRDS(pbmc, file = "CellType_pbmc.rds")
pbmc@meta.data$celltype_HPCA = "NA"
for(i in 1:nrow(celltype)){
  pbmc@meta.data[which(pbmc@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype_HPCA'] <- celltype$celltype[i]}
p1 = DimPlot(pbmc, group.by="celltype_HPCA", repel=T, label=T, label.size=5, reduction='tsne')
p2 = DimPlot(pbmc, group.by="celltype_HPCA", repel=T, label=T, label.size=5, reduction='umap')
p3 = p1+p2+ plot_layout(guides = 'collect')
ggsave("tSNE_celltype_HPCA.png", p1, width=7 ,height=6)
ggsave("UMAP_celltype_HPCA.png", p2, width=7 ,height=6)
ggsave("celltype_HPCA.png", p3, width=10 ,height=5)
#====================================================================================================

#人工注释
#Alpha
c1=VlnPlot(pbmc, features = c("PLCE1"))
ggsave("Alpha1.png", c1, width=5 ,height=5)
#Beta
c2=VlnPlot(pbmc, features = c("INS"))
ggsave("Beta1.png", c2, width=5 ,height=5)
#Delta
c3=VlnPlot(pbmc, features = c("PDLIM4"))
ggsave("Delta1.png", c3, width=5 ,height=5)
#Gamma\pp
c4=VlnPlot(pbmc, features = c("PPY"))
ggsave("Gamma1.png", c4, width=5 ,height=5)
#Acinar
c5=VlnPlot(pbmc, features = c("PRSS1"))
ggsave("Acinar1.png", c5, width=5 ,height=5)
#Ductal
c6=VlnPlot(pbmc, features = c("CFTR"))
ggsave("Ductal1.png", c6, width=5 ,height=5)
#Stellate
#c7=VlnPlot(pbmc, features = c("COL3A1", "BGN", "IGFBP5"))
#ggsave("Stellate1.png", c7, width=10 ,height=6)



#Alpha
c11=FeaturePlot(pbmc, features = c("PLCE1"),reduction = "tsne")#tsne的映射图
#c11=FeaturePlot(pbmc, features = c("PLCE1"))#umap的映射图
ggsave("Alpha22.png", c11, width=5 ,height=5)

#Beta
c22=FeaturePlot(pbmc, features = c("INS"),reduction = "tsne")
ggsave("Beta22.png", c22, width=5 ,height=5)
#Delta
c33=FeaturePlot(pbmc, features = c("PDLIM4"),reduction = "tsne")
ggsave("Delta22.png", c33, width=5 ,height=5)
#Gamma
c44=FeaturePlot(pbmc, features = c("PPY"),reduction = "tsne")
ggsave("Gamma22.png", c44, width=5 ,height=5)
#Acinar
c55=FeaturePlot(pbmc, features = c("PRSS1"),reduction = "tsne")
ggsave("Acinar22.png", c55, width=5 ,height=5)
#Ductal
c66=FeaturePlot(pbmc, features = c("CFTR"),reduction = "tsne")
ggsave("Ductal22.png", c66, width=5 ,height=5)
#Stellate
#c77=FeaturePlot(pbmc, features = c("COL3A1", "BGN", "IGFBP5"))
#ggsave("Stellate2.png", c77, width=10 ,height=6)


## =============cell type（人工注释）
new.cluster.ids <- c("β细胞",   
                     "α细胞",    
                     "腺管细胞",   
                     "PP细胞",           
                     "腺泡细胞",         
                     "Unidentified",  
                     "δ细胞",             
                     "Unidentified")  
names(new.cluster.ids) <- levels(pbmc)
new.cluster.ids
# before
Idents(pbmc)
# rename
pbmc <- RenameIdents(pbmc, new.cluster.ids)
Idents(pbmc)
head(pbmc@meta.data)

#可视化
p3=DimPlot(pbmc, reduction = "umap", label = TRUE, label.size=5,pt.size = 1.2) + NoLegend()
p4=DimPlot(pbmc, reduction = "tsne", label = TRUE,label.size=5, pt.size = 1.2) + NoLegend()
ggsave("人工_umap.png",p3,width = 10,height = 8)
ggsave("人工_tsne.png",p4,width = 10,height = 8)

# save
pbmc@meta.data$cell_anno <- Idents(pbmc)
write.csv(pbmc@meta.data,file = "metadata.csv")
saveRDS(pbmc, file = "pbmc3k_final.rds")
newcelltype = data.frame(names(new.cluster.ids),new.cluster.ids, stringsAsFactors = F)
pbmc@meta.data$newcelltype = "NA"
for(i in 1:nrow(newcelltype)){
  pbmc@meta.data[which(pbmc@meta.data$seurat_clusters == newcelltype$names.new.cluster.ids.[i]),'newcelltype'] <- newcelltype$new.cluster.ids[i]}
p1 = DimPlot(pbmc, group.by="newcelltype", label=T, label.size=5, reduction='tsne')
p2 = DimPlot(pbmc, group.by="newcelltype", label=T, label.size=5, reduction='umap')
p3 = plotc <- p1+p2+ plot_layout(guides = 'collect')
ggsave("人工注释.tSNE_celltype.png", p1, width=10 ,height=7)
ggsave("人工注释.UMAP_celltype.png", p2, width=10 ,height=7)
ggsave("人工注释.celltype.pdf", p3, width=15 ,height=7)
ggsave("人工注释.celltype.png", p3, width=15 ,height=6)
