%GKTODO: correct path

atlaspath = 'C:\Users\hgaddy1\MATLAB\MNI152_T1_1mm_brain_labels_cropped.nii';
datapath = 'C:\Users\hgaddy1\MATLAB\ADHD200_40sub_preprocessed\ForGreg\';

%% Load data - some hardcoding in this prototype

% load ROIs
nn = load_nii(atlaspath);
roiLabel = nn.img;
roiLabel(roiLabel > 100) = roiLabel(roiLabel > 100) - 65;
roiLabelMask = roiLabel;

cd(datapath);
files = dir('*.nii.gz');
%% Compute matrix
% matrix is 40x70: each row represents sum, of average BOLD levels for
% each of the 70 regions. 
% Equivalently, each column represents average BOLD levels for each 
% participant.

vecs = zeros(length(files), 70);
for f=1:length(files)
    %read it
    temp_img = load_nii(files(f),name);
    fmri = temp_img.img;
    
    %average over time
    avg = mean(fmri,4);
    
    %assign each voxel intensity to a label
    for i=1:size(avg,1)
        for j=1:size(avg,2)
           for k=1:size(avg,3) 
                region = roiLabelMask(i,j,k);
                vecs(f,region) = vecs(f,region) + avg(i,j,k);
           end
        end
    end 
    
end


save('partroisum_matrix.mat',vecs)
%now you have 1x70 for each person

%% if wanting to do more with time series: consider adding the following

% temp = importdata('kki42_subjectinformation.csv');
% 
% for i = 2:length(temp)
%     reorderIdx(i-1) = str2num(temp{i}(end-1:end));
% end

% for t=1:size(fmri, 4)
%     f = fmri(:,:,:,t);
%     
% end
