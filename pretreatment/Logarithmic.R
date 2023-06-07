rm(list = ls())
dat<-read.csv("E:\\毕业论文\\单细胞数据集\\GSE83139\\GSE83139（降序后）.csv",header = F)
dat <- dat[-1,] # 去掉第一行
dat <- dat[,-1]# 去掉列名
dat2=as.data.frame(lapply(dat,as.numeric))#数据框数据转换为数值型
ndat <- log(1+dat2)
write.csv(ndat,file ="GSE83139（降序且取对数）.csv")
