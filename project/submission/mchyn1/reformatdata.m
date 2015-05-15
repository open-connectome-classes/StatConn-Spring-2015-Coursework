clear all;close all;clc;
tic

%% Paths

atlaspath = '/MNI152_T1_1mm_brain_labels_cropped.nii';

%% Loading data

map=load_nii(atlaspath);
roi=map.img;
roi(roi>100)=roi(roi>100)-65;

cd('New Folder');
files=dir('*.nii.gz');
cd ..;
info=csvread('ADHDlabels_ordered.csv',1,0);
info(:,4)=1:40;

%% Formatting Data
% ends up with ROI x T size matrix stored in variable m
% rlabel is assigned with the same index to keep track of region

m=cell(40,1);
M=cell(21,1);
for i = 5 %length(files)
    S = load_nii(files(i).name);
    data = S.img;
    dims=size(data);
    m{i}=zeros(70,dims(4));
    region=nan(1);
    count=zeros(70,1);
    temp=zeros(1,dims(4));
    rlabel=nan(70,1);
    for x = 1:dims(1)-1
        for y = 1:dims(2)-1
            for z = 1:dims(3)-1
                if roi(x,y,z)~=0
                    region=roi(x,y,z);
                    count(region)=count(region)+1;
                    for t=1:dims(4)
                        m{i}(region,t)=m{i}(region,t)+data(x,y,z,t);
                    end
                end
            end
        end
    end
    m{i}(region,t)=m{i}(region,t)./count(region);
    hasdata=find(count); %contains regions with data
    for j=1:length(hasdata)
        M{i}(j,:)=m{i}(hasdata(j),:);
    end
end
toc
%save('regbytime',M);
