clear
data=xlsread('E:\毕业论文\单细胞数据集\GSE83139\网络构建\GSE83139beta-child.xlsx');
alpha=0.1;
boxsize=1.5;
weighted=1;
c = [];
kk=1;
cndm = condition_ndm(data,alpha,boxsize,kk);
xlswrite('E:\毕业论文\单细胞数据集\GSE83139\网络构建\度矩阵\GSE83139beta-child-cndm.xlsx',cndm);