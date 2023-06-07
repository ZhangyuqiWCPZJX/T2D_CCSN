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

[n1,n2] = size(data);%n1=��,n2=��
upbound = zeros(n1,n2);%���ϵ�
lowbound = zeros(n1,n2);%����

for i = 1 : n1
    [s1,s2] = sort(data(i,:));%ȡ��i�н�������
    n3 = n2-sum(sign(s1));%δ���ϸ��
    h = round(boxsize*sqrt(sum(sign(s1)))); %��ֵ����һ�����е�ĸ���
    k = 1;
    while k <= n2
        s = 0;
        while k+s+1 <= n2 && s1(k+s+1) == s1(k)%s�Ƿ��ʾ�������ͬ��ϸ����������
            s = s+1;
        end
        if s >= h
            upbound(i,s2(k:k+s)) = data(i,s2(k));%����
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
    adjmc = (c*n2-a*a')./sqrt((a*a').*((n2-a)*(n2-a)')/(n2-1)+eps);  %�������4s-3
    adjmc = (adjmc > p);
    if kk ~=0 %�����ڣ�����1
    id = condition_g(adjmc,kk);
    adjmc = zeros(n1,n1);%������һ��n1��n1�е�ȫ�����
    for m = 1:kk   
    B_z = bsxfun(@times,B(id(m),:),B);
    idc = find(B(id(m),:)~=0);%���ؾ����з���Ԫ�ص�λ��
    B_z = B_z(:,idc);
    [~,r] = size(B_z);%�����������
    a_z = sum(B_z,2);
    c_z = B_z*B_z';
    adjmc1 =(c_z*r-a_z*a_z')./sqrt((a_z*a_z').*((r-a_z)*(r-a_z)')/(r-1)+eps);
    adjmc1 = (adjmc1 > p);%p=2.3263
    adjmc = adjmc + adjmc1;%��Ӧλ�����
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