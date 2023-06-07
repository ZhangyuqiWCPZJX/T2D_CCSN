rm(list=ls())
setwd("C:\\Users\\Y406\\Documents\\暗基因")
library(ggplot2)
theme_set(theme_bw())#为白色背景主题
# theme_grey()为默认主题,theme_bw()为白色背景主题,theme_classic()为经典主题
df1<-read.csv("E:\\毕业论文\\单细胞数据集\\GSE83139\\暗基因\\暗基因线图\\MRAS\\MRAS基因.csv",header = T)
df1$stage<-factor(df1$stage,levels = c("health","T2D"))
df2<-read.csv("E:\\毕业论文\\单细胞数据集\\GSE83139\\暗基因\\暗基因线图\\MRAS\\MRAS度矩阵.csv",header = T)
df2$stage<-factor(df2$stage,levels = c("health","T2D"))
df11<-read.csv("E:\\毕业论文\\单细胞数据集\\GSE83139\\暗基因\\暗基因线图\\MRAS\\MRAS基因平均.csv",header = T)
df11$stage<-factor(df11$stage,levels = c("health","T2D"))
df22<-read.csv("E:\\毕业论文\\单细胞数据集\\GSE83139\\暗基因\\暗基因线图\\MRAS\\MRAS度平均.csv",header = T)
df22$stage<-factor(df22$stage,levels = c("health","T2D"))

#散点加折线
g<-ggplot(df1, aes(x=stage, y=GEM, color = stage))
g+geom_point(position = "jitter", alpha=.3)+
  geom_line(data = df11,aes(x=stage, y=average),color="black",group=df11$group,size=0.5) +
  ylim(c(0,2))+
  labs(y="EXP", 
       x="stage", 
       title="MRAS")
g<-ggplot(df2, aes(x=stage, y=NDM, color = stage))
g+geom_point(position = "jitter", alpha=.3)+
  geom_line(data = df22,aes(x=stage, y=average),color="black",group=df22$group,size=0.5) +
  ylim(c(0,2))+
  labs(y="Degree", 
       x="stage", 
       title="MRAS")

#折线
ggplot(df11, aes(stage,average,group=group)) +
  geom_point() +
  geom_line() +
  labs(y="EXP", 
       x="stage", 
       title="XIST")
ggplot(df22, aes(stage,average,group=group)) +
  geom_point() +
  geom_line() +
  labs(y="Degree", 
       x="stage", 
       title="XIST")

#散点
g<-ggplot(df1, aes(x=stage, y=GEM, color = stage))
g+geom_point(position = "jitter", alpha=.3)+
  labs(y="EXP", 
       x="stage", 
       title="XIST")
g<-ggplot(df2, aes(x=stage, y=GEM, color = stage))
g+geom_point(position = "jitter", alpha=.3)+
  labs(y="Degree", 
       x="stage", 
       title="XIST")
