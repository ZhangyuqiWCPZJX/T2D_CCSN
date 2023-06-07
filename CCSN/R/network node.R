if(!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('WGCNA')
library(WGCNA)

TOM <-read.csv('E:\\毕业论文\\单细胞数据集\\GSE83139\\网络构建\\网络\\GSE83139alpha-T2D-CCSN(avr).csv',header=T,row.names = 1)

cyt1= exportNetworkToCytoscape(TOM,threshold = 0.98,
                               edgeFile="E:\\毕业论文\\单细胞数据集\\GSE83139\\网络构建\\网络\\网络节点\\alpha-T2D\\edge0.98_alpha-t2d.txt",
                               nodeFile="E:\\毕业论文\\单细胞数据集\\GSE83139\\网络构建\\网络\\网络节点\\alpha-T2D\\node0.98_alpha-t2d.txt",weighted = TRUE)


#扩内存
memory.size(F)
memory.limit()
memory.limit(size=10000)