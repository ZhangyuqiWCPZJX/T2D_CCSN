function ccsn = csnet(data,c,alpha,boxsize,kk,weighted)
%Define the neighborhood of each plot
if nargin < 6 || isempty(weighted)
    weighted = 0;
end
if nargin < 4 || isempty(boxsize)
    boxsize = 0.1;
end
if nargin <3 || isempty(alpha)
    alpha = 0.01;
end
 
[n1,n2] = size(data);
if nargin <2 || isempty(c)
    c = 1 : n2;
end

[n1,n2] = size(data);%n1=行,n2=列
upbound = zeros(n1,n2);%向上的
lowbound = zeros(n1,n2);%下限

for i = 1 : n1
    [s1,s2] = sort(data(i,:));%取第i行进行升序
    n3 = n2-sum(sign(s1));%未表达细胞
    h = round(boxsize*sqrt(sum(sign(s1)))); %阈值，第一个框中点的个数
    k = 1;
    while k <= n2
        s = 0;
        while k+s+1 <= n2 && s1(k+s+1) == s1(k)%s是否表示表达量相同的细胞个数？？
            s = s+1;
        end
        if s >= h
            upbound(i,s2(k:k+s)) = data(i,s2(k));%不懂
            lowbound(i,s2(k:k+s)) = data(i,s2(k));
        else
            upbound(i,s2(k:k+s)) = data(i,s2(min(n2,k+s+h)));
            lowbound(i,s2(k:k+s)) = data(i,s2(max(n3*(n3>h)+1,k-h)));
        end
        k = k+s+1;
    end
end
 
%Construction of cell-specific network
ccsn = cell(1,n2);
B = zeros(n1,n2);
p = -icdf('norm',alpha,0,1);
for k = 1:n2
    for j = 1 : n2
       B(:,j) = (data(:,j) <= upbound(:,k) & data(:,j) >= lowbound(:,k)) & data(:,k);
    end
    a = sum(B,2);
    c = B*B';
    adjmc = (c*n2-a*a')./sqrt((a*a').*((n2-a)*(n2-a)')/(n2-1)+eps);  %补充材料4s-3
    adjmc = (adjmc > p);
    if kk ~=0 %不等于，返回1
    id = condition_g(adjmc,kk);
    adjmc = zeros(n1,n1);%生成了一个n1行n1列的全零矩阵
    for m = 1:kk   
    B_z = bsxfun(@times,B(id(m),:),B);
    idc = find(B(id(m),:)~=0);%返回矩阵中非零元素的位置
    B_z = B_z(:,idc);
    [~,r] = size(B_z);%忽略输出参数
    a_z = sum(B_z,2);
    c_z = B_z*B_z';
    adjmc1 =(c_z*r-a_z*a_z')./sqrt((a_z*a_z').*((r-a_z)*(r-a_z)')/(r-1)+eps);
    adjmc1 = (adjmc1 > p);%p=2.3263
    adjmc = adjmc + adjmc1;%对应位置相加
    end
    else 
        kk=1;
    end
    
    if weighted
        ccsn{k} = adjmc.*(adjmc > 0);
    else
        ccsn{k} = sparse(adjmc > p);
    end
    disp(['Cell ' num2str(k) ' is completed']);
end