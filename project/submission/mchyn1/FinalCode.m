%import data
url='https://www.dropbox.com/sh/idt3d0gylplyo31/AADtmR1GhccrhCG7RY6C1O5Qa/mchyn1_aceeccc.zip?dl=1';
cmd=['wget ' url ' -O "mchyn1_aceeccc.zip" --no-check-certificate'];
system(cmd, '-echo');
unzip('mchyn1_aceeccc.zip');
addpath('mchyn1_aceeccc');

%not caring about the pairs of age/gender, sorting by phenotype
in1 = 1;
in2 = 1;
adhd=nan(20,1);
tdc=nan(20,1);
for i = 1:40
    if info(i,3) == 1
        adhd(in1)=info(i,4);
        in1=in1+1;
    else
        tdc(in2)=info(i,4);
        in2=in2+1;
    end
    
end
%ADHD
A1=cell(20,1);
for i = 1:20
    A1{i}=A{adhd(i)};
end
A1 = A1(~cellfun(@isempty, A1));

A2=cell(20,1);
for i =1:20
    A2{i}=A{tdc(i)};
end
A2=A2(~cellfun(@isempty,A2));

%kernel
K=nan(21,1);
X=cell(21,1);
%(21,14,21); %voxel, patient, v1v1 v1v2 v1v3 etc
Anew=[A1;A2];
figure
for i = 1:length(Anew)
    subplot(4,4,i);imagesc(Anew{i});
    colorbar
end
suptitle('region by region autoregression matrices. top 2 rows ADHD, bottom 2 TDC');
%TDC
for i = 1:14
    for v = 1:21
        X{v}(i,:)=Anew{i}(v,:); %X{region}(x1,x2,x3..) where x1=v1v1 v1v2..
    end
end
for v = 1:21
    K(v)=exp(-norm(X{v}(1:7,:)-X{v}(8:14,:))).^2;
end
figure
imagesc(X{4});
colorbar
xlabel('region of interest');
ylabel('patient');
title('4th row of A for each patient');
%permutation testing
Kperm=nan(21,10000);
for n=1:10000
    p = randperm(14); %permuting labels
    Aperm=cell(14,1);
    Xperm=cell(21,1);
    for i = 1:14
        Aperm{i}=Anew{p(i)};
    end
    for i = 1:14
        for v = 1:21
            Xperm{v}(i,:)=Aperm{i}(v,:);
        end
    end
    for v = 1:21
        Kperm(v,n)=exp(-norm(Xperm{v}(1:7,:)-Xperm{v}(8:14,:))).^2;
    end
end

%significance using .05 level
count=zeros(21,1);
for n = 1:10000
    for i = 1:21
        if Kperm(i,n)>K(i)
            count(i)=count(i)+1;
        end
    end
end
count=1-count/10000;
figure
plot(1:21,count,'.');
hold on
sig=nan(21,1);
for i = 1:21
    if count(i)<.05
        sig(i)=1;
        plot(i,count(i),'*r');
    else
        sig(i)=0;
    end
end
xlabel('region of interest');
ylabel('p-value');
title('.05 level significance of region');