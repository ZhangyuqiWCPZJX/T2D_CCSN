install.packages('BiocManager') #需要首先安装 BiocManager，如果尚未安装请先执行该步
BiocManager::install('edgeR')

library(edgeR)
#setwd("C:\\Users\\Y406\\Documents")
exprSet_all <- read.csv("E:\\毕业论文\\单细胞数据集\\GSE83139\\网络构建\\基因表达矩阵\\beta\\GSE83139beta-T2D-child.csv")
exprSet <- exprSet_all[,-1]
group_list <- factor(c(rep("T2D.beta",38),rep("child.beta",19)))
exprSet <- DGEList(counts = exprSet, group = group_list)
exprSet <- calcNormFactors(exprSet)
exprSet <- estimateCommonDisp(exprSet)
exprSet <- estimateTagwiseDisp(exprSet)
et <- exactTest(exprSet)
tTag <- topTags(et, n=nrow(exprSet))
tTag <- as.data.frame(tTag)
write.csv(tTag,file = "E:\\毕业论文\\单细胞数据集\\GSE83139\\差异基因\\T2Dvschild.beta.csv")
