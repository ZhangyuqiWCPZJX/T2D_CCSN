clear
data=xlsread('E:\��ҵ����\��ϸ�����ݼ�\GSE83139\���繹��\GSE83139beta\GSE83139beta-T2D.xlsx');
alpha=0.1;
boxsize=1.5;
weighted=1;
c = [];
kk=1;
ccsn= csnet(data,c,alpha,boxsize,kk,weighted);

data = zeros(14999);
for i = 1:length(ccsn)
    tmp = ccsn{1,i};
    data = data + tmp;    
end
avr = data/length(ccsn);

csvwrite('E:\��ҵ����\��ϸ�����ݼ�\GSE83139\���繹��\GSE83139beta\GSE83139beta-T2D-CCSN(avr).csv',avr);