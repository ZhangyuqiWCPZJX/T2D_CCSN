clear
data=xlsread('E:\��ҵ����\��ϸ�����ݼ�\GSE83139\���繹��\GSE83139beta-child.xlsx');
alpha=0.1;
boxsize=1.5;
weighted=1;
c = [];
kk=1;
cndm = condition_ndm(data,alpha,boxsize,kk);
xlswrite('E:\��ҵ����\��ϸ�����ݼ�\GSE83139\���繹��\�Ⱦ���\GSE83139beta-child-cndm.xlsx',cndm);