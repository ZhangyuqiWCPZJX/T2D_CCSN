rm(list=ls())
setwd("E:\\毕业论文\\单细胞数据集\\GSE83139\\富集分析")
install.packages('GOplot')
install.packages('ggnewscale')
BiocManager::install('topGO')
BiocManager::install('ComplexHeatmap')
library(openxlsx)#读取.xlsx文件
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例

#读取差异表达基因，将基因ID从GENE_SYMBOL转换为ENTREZ_ID：
#载入差异表达数据，只需基因ID(GO,KEGG,GSEA需要)和Log2FoldChange(GSEA需要)即可
info <- read.xlsx( "E:\\毕业论文\\单细胞数据集\\GSE83139\\暗基因\\alpha-childvsT2D暗.xlsx", rowNames = F,colNames = T)

#指定富集分析的物种库
GO_database <- 'org.Hs.eg.db' 
#GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
KEGG_database <- 'hsa' 
#KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html

#gene ID转换
gene <- bitr(info$gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)

#GO分析
erich.go.BP = enrichGO(gene = gene$ENTREZID,
                       OrgDb = GO_database,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
##分析完成后，作图
dotplot(erich.go.BP)
barplot(erich.go.BP)#画柱形图
#plotGOgraph(erich.go.BP)树形图

erich.go.CC = enrichGO(gene = gene$ENTREZID,
                       OrgDb = GO_database,
                       keyType = "ENTREZID",
                       ont = "CC",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
##分析完成后，作图
dotplot(erich.go.CC)
barplot(erich.go.CC)
# plotGOgraph(erich.go.cc)

erich.go.MF = enrichGO(gene = gene$ENTREZID,
                       OrgDb = GO_database,
                       keyType = "ENTREZID",
                       ont = "MF",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
##分析完成后，作图
dotplot(erich.go.MF)
barplot(erich.go.MF)
# plotGOgraph(erich.go.MF)


#KEGG分析
KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.2,
                 qvalueCutoff = 0.2)
dotplot(KEGG)
barplot(KEGG)
# plotGOgraph(KEGG)


hsa_kegg<-download_KEGG("hsa")
str(hsa_kegg)
length(unique(hsa_kegg$KEGGPATHID2EXTID$from)) #信号通路个数345

length(unique(hsa_kegg$KEGGPATHID2EXTID$to))  #所有信号通路中的基因数8112

x<-data.frame(hsa_kegg$KEGGPATHID2EXTID)  #两列为pathwayid,基因id extid 扩增标识符
#一个信号通路含有多个基因，一个基因在多个信号通路

length(unique(hsa_kegg$KEGGPATHID2NAME$from)) ###信号通路对应的名称548
length(unique(hsa_kegg$KEGGPATHID2NAME$to)) ###信号通路对应的名字548
y<-data.frame(hsa_kegg$KEGGPATHID2NAME)  ###信号通路信息 

info<-gene
#ID转换
deg.id <- bitr(info$SYMBOL, #基因名
               fromType = "SYMBOL", #从gene symbol
               toType = "ENTREZID", #提取ENTREZ ID
               OrgDb = "org.Hs.eg.db") #相应物种的包
#合并基因名、ENTREZID
idvec <- deg.id$ENTREZID
names(idvec) <- deg.id$SYMBOL
info$ENTREZID <- idvec[info$SYMBOL]
#保存到文件
write.csv(info, "alpha(child_t2d)暗基因ENTREZID.csv", quote = F, row.names = F)
#deg <- read.csv("diff_ENTREZID.csv", as.is = T)


#KEGG富集
ekk <- enrichKEGG(gene = info$ENTREZID,
                  keyType = 'kegg',
                  organism = 'hsa',
                  pvalueCutoff = 0.2,
                  pAdjustMethod  = "BH"
                  #qvalueCutoff  = 0.05
)
dim(ekk)
#8 9
#把富集分析结果保存到文件
write.csv(ekk,"enrichKEGG.csv",quote = F)
#可以用语义学方法，合并相似的KEGG term。需要较长时间。
#参考资料：<https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/>
#ekk_1 <- simplify(ekk, cutoff = 0.7, by = "p.adjust", select_fun = min)
#dim(ekk_1)
#write.csv(ekk_1,"enrichKEGG_simple.csv",quote = F)

dotplot(ekk,showCategory=10)
barplot(ekk)
#对于基因和富集的pathways之间的对应关系进行展示，如果一个基因位于一个pathway下，则将该基因与pathway连线
cnetplot(ekk, categorySize="pvalue", showCategory = 5)

###upsetplot
#BiocManager::install("ggupset")
#BiocManager::install("pathview")
library(ggupset)
upsetplot(ekk)

library(pathview)
#Downloading xml/png files for hsa00190 会比较慢
pv.out <- pathview(gene.data = info$ENTREZID,#需要提供的基因向量，默认是Entrez_ID。可以由gene.idtype设置
                   pathway.id = 'hsa00190', #在KEGG中的ID
                   species ="hsa",
                   #out.suffix = "gse57691"，#文件输出名称，默认是pathview
                   kegg.native =TRUE#默认是TRUE输出完整pathway的png格式文件，反之输出pdf文件
)
#在pathway通路图上标记富集到的基因
browseKEGG(ekk, "hsa00190")
#https://www.kegg.jp/kegg-bin/show_pathway?hsa00190/27089/55967/522/4706/4716/1345/4702/539/4694/1346/1350/1329/4725/1347/10063/1340


# #GSEA分析
# names(info) <- c('SYMBOL','Log2FoldChange','pvalue','padj')
# info_merge <- merge(info,gene,by='SYMBOL')#合并转换后的基因ID和Log2FoldChange
# GSEA_input <- info_merge$Log2FoldChange
# names(GSEA_input) = info_merge$ENTREZID
# GSEA_input = sort(GSEA_input, decreasing = TRUE)
# GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = 0)#GSEA富集分析              


## GO/KEGG富集柱状图+点状图(有问题):
GO<-enrichGO( gene$ENTREZID,#GO富集分析
              OrgDb = GO_database,
              keyType = "ENTREZID",#设定读取的gene ID类型
              ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
              pvalueCutoff = 0.05,#设定p值阈值
              qvalueCutoff = 0.05,#设定q值阈值
              readable = T)
KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.2,
                 qvalueCutoff = 0.2)
barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#柱状图
barplot(KEGG,showCategory = 40,title = 'KEGG Pathway')
dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#点状图
dotplot(KEGG)


# 富集基因与所在功能集/通路集的关联网络图
enrichplot::cnetplot(GO,circular=TRUE,colorEdge = TRUE)#基因-通路关联网络图
enrichplot::cnetplot(KEGG,circular=TRUE,colorEdge = TRUE)#circluar为指定是否环化，基因过多时建议设置为FALSE

# 也可以以热图形式展现关联关系:
enrichplot::heatplot(GO,showCategory = 50)#基因-通路关联热图
enrichplot::heatplot(KEGG,showCategory = 50)

#富集到的功能集/通路集之间的关联网络图
GO2 <- pairwise_termsim(GO)
KEGG2 <- pairwise_termsim(KEGG)
enrichplot::emapplot(GO2,showCategory = 50, color = "p.adjust", layout = "kk")#通路间关联网络图
enrichplot::emapplot(KEGG2,showCategory =50, color = "p.adjust", layout = "kk")


#保存KEGG富集到的通路至本地文件并选择通路进行展示
write.table(KEGG$ID, file = "KEGG_IDs.txt", #将所有KEGG富集到的通路写入本地文件查看
            append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
browseKEGG(KEGG,"hsa04666")#选择其中的hsa04666通路进行展示


#GO富集功能网络图
GO_BP<-enrichGO( gene$ENTREZID,#GO富集分析BP模块
                 OrgDb = GO_database,
                 keyType = "ENTREZID",
                 ont = "BP",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 minGSSize = 10,
                 maxGSSize = 500,
                 readable = T)
plotGOgraph(GO_BP)#GO-BP功能网络图
GO_CC<-enrichGO( gene$ENTREZID,#GO富集分析CC模块
                 OrgDb = GO_database,
                 keyType = "ENTREZID",
                 ont = "CC",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 minGSSize = 10,
                 maxGSSize = 500,
                 readable = T)
plotGOgraph(GO_CC)#GO-CC功能网络图
GO_MF<-enrichGO( gene$ENTREZID,#GO富集分析MF模块
                 OrgDb = GO_database,
                 keyType = "ENTREZID",
                 ont = "MF",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 minGSSize = 10,
                 maxGSSize = 500,
                 readable = T)
plotGOgraph(GO_MF)#GO-MF功能网络图




#GSEA分析
info <- read.xlsx( "C:\\Users\\Y406\\Documents\\内分泌细胞暗基因.xlsx", rowNames = F,colNames = T)
names(info) <- c('SYMBOL','Log2FoldChange','PValue','FDR')
info_merge <- merge(info,gene,by='SYMBOL')#合并转换后的基因ID和Log2FoldChange
GSEA_input <- info_merge$Log2FoldChange
names(GSEA_input) = info_merge$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)
GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = 0)#GSEA富集分析


#GO富集弦图
genedata<-data.frame(ID=info$SYMBOL,logFC=info$Log2FoldChange)
write.table(GO$ONTOLOGY, file = "GO_ONTOLOGYs.txt", #将所有GO富集到的基因集所对应的类型写入本地文件从而得到BP/CC/MF各自的起始位置如我的数据里是1，2103，2410
            append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
GOplotIn_BP<-GO[1:10,c(2,3,7,9)] #提取GO富集BP的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_CC<-GO[1:14,c(2,3,7,9)]#提取GO富集CC的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_MF<-GO[619:629,c(2,3,7,9)]#提取GO富集MF的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_BP$geneID <-str_replace_all(GOplotIn_BP$geneID,'/',',') #把GeneID列中的’/’替换成‘,’
GOplotIn_CC$geneID <-str_replace_all(GOplotIn_CC$geneID,'/',',')
GOplotIn_MF$geneID <-str_replace_all(GOplotIn_MF$geneID,'/',',')
names(GOplotIn_BP)<-c('ID','Term','adj_pval','Genes')#修改列名,后面弦图绘制的时候需要这样的格式
names(GOplotIn_CC)<-c('ID','Term','adj_pval','Genes')
names(GOplotIn_MF)<-c('ID','Term','adj_pval','Genes')
GOplotIn_BP$Category = "BP"#分类信息
GOplotIn_CC$Category = "CC"
GOplotIn_MF$Category = "MF"
circ_BP<-GOplot::circle_dat(GOplotIn_BP,genedata) #GOplot导入数据格式整理
circ_CC<-GOplot::circle_dat(GOplotIn_CC,genedata) 
circ_MF<-GOplot::circle_dat(GOplotIn_MF,genedata) 
chord_BP<-chord_dat(data = circ_BP,genes = genedata) #生成含有选定基因的数据框
chord_CC<-chord_dat(data = circ_CC,genes = genedata) 
chord_MF<-chord_dat(data = circ_MF,genes = genedata) 
GOChord(data = chord_BP,#弦图
        title = 'GO-Biological Process',space = 0.01,#GO Term间距
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), #上下调基因颜色
        process.label = 10) #GO Term字体大小
GOChord(data = chord_CC,title = 'GO-Cellular Component',space = 0.01,
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), 
        process.label = 10) 
GOChord(data = chord_MF,title = 'GO-Mollecular Function',space = 0.01,
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), 
        process.label = 10)



#GO富集弦表图
GOCircle(circ_BP) #弦表图
GOCircle(circ_CC) 
GOCircle(circ_MF) 


#GO富集系统聚类图
chord<-chord_dat(data = circ_BP,genes = genedata) #生成含有选定基因的数据框
GOCluster(circ_BP,GOplotIn_BP$Term) #系统聚类图
chord<-chord_dat(data = circ_CC,genes = genedata)
GOCluster(circ_CC,GOplotIn_CC$Term) 
chord<-chord_dat(data = circ_MF,genes = genedata) 
GOCluster(circ_MF,GOplotIn_MF$Term) 



# #GSEA富集图
# ridgeplot(GSEA_KEGG) 
# gseaplot2(GSEA_KEGG,1)
# gseaplot2(GSEA_KEGG,1:30)#30是根据ridgeplot中有30个富集通路得到的



