%% eval_partitions.m
% MRI Partition Comparion for small graphs generated on the KKI2009,
% 21-subject 42-scan dataset. 
%

function performance = eval_partitions(metric, N)

%Sets up number of subjects and similarity metric
if ~exist('N', 'var')
    N = 42;
end
if ~exist('metric', 'var') %modes of the 'norm' function, for now
    metric = 'fro';
    %     metric = 2;
    %     metric = 1;
end

%gets directory structure for files to read in
dirinfo = dir();
dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
subdirinfo = cell(length(dirinfo));
for K = 1 : length(dirinfo)
    thisdir = dirinfo(K).name;
    subdirinfo{K} = dir(fullfile(thisdir, '*.mat'));
end

%reorders subject based on indx
temp = importdata('kki42_subjectinformation.csv');
for i = 2:length(temp)
    reorderIdx(i-1) = str2num(temp{i}(end-1:end));
end

cwd = pwd;

w=0;
figure(2);
%does TRT comparison for each dataset
for ii=1:length(subdirinfo)
    if length(subdirinfo{ii}) ~= N
        continue; %ignores . and .. directories
    end
    dirinfo(ii).name
    cd(dirinfo(ii).name); %move to subdirectory
    w = w+1;
    
    %load graphs in subject order
    c = 1;
    for i = reorderIdx
        temp = load(subdirinfo{ii}(i).name);
        tgraph = log10(full(temp.graph));
        tgraph(isinf(tgraph))=0;
        sg(:,:,c) = tgraph;
        c = c+1;
    end
    
    %compute similarity
    for i = 1:size(sg,3)
        for j = 1:size(sg,3)
            gErr(i,j) = norm(sg(:,:,i)-sg(:,:,j), metric);
        end
    end
    
    %compute TRT
    matches = 0;
    for i = 1:size(gErr,1)
        temp = sort(gErr(i,:));
        q = i-1+2*mod(i,2);
        matches = matches + (temp(2)==gErr(i,q));
    end
    
%     strcat(num2str(matches), ' / ', num2str(N))

    %save result
    performance(w).atlas = dirinfo(ii).name;
    performance(w).matches = matches;
    performance(w).gErr = gErr;
    
    %show result
    subplot(2,2,w);
    imagesc(gErr); caxis([0, 200]);
    title(strcat(dirinfo(ii).name, ' (N=', num2str(size(sg,1)), ') atlas: ', num2str(matches), '/',num2str(N)));
        
    %reset to baseline
    cd(cwd);
    clear sg
    clear gErr
    clear matches
end

end