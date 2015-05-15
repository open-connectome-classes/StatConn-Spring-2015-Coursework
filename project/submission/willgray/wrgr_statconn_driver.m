%% Stat Conn
% W. Gray Roncal

% Get data
!wget https://www.dropbox.com/sh/idt3d0gylplyo31/AAAFirAczHNiqLar4UHkXJkPa/wrgr_graphs.zip
unzip('wrgr_graphs.zip')
%% Download from OCP
oo = OCP();
oo.setServerLocation('braingraph1dev.cs.jhu.edu');
oo.setAnnoToken('img2graph_v29_gala_node_paramset_5');
oo.setImageToken('kasthuri11cc');
q = OCPQuery;
q.setType(eOCPQueryType.annoDense);
q.setCutoutArgs([5472,6496],[8712,9736],[1000,1100],1);
annoDetect = oo.query(q);

q.setType(eOCPQueryType.imageDense);
imData = oo.query(q);

q.setType(eOCPQueryType.annoDense);
oo.setAnnoToken('ac3');
annoTruth = oo.query(q);

oo.setServerLocation('openconnecto.me');
oo.setAnnoToken('ac3_synTruth_v4');
synTruth = oo.query(q);

oo.setServerLocation('braingraph1dev.cs.jhu.edu');
oo.setAnnoToken('img2graph_v29_syn_obj_edge_paramset_10');
synTest = oo.query(q);

h = image(imData); h.associate(synTest)

%% Compute line graphs and error

% Preproc
temp = cleanup_speckle(annoTruth.data);
annoTruth.setCutout(temp);

% Associate Synapses
eTruth = statconn_associate_neusyn(annoTruth, synTruth);
eTest = statconn_associate_neusyn(annoDetect, synTruth);

% Build Graphs
[~,~, lgTruth, ~] = statconn_construct_graph(eTruth);
[~,~, lgAnnoDetect, ~] = statconn_construct_graph(eTest);

% Compute Error
gErrAnnoDetect = statconn_compute_graph_error_simple(lgAnnoDetect, lgTruth);

sprintf('F1 Score for best Gala operating point: %0.2f\n',gErrAnnoDetect.scoreF1)
%% Find orphans

% Find things that hit boundaries
anno = annoDetect.data;
uannoEdge = unique([unique(anno(:,:,1)); unique(anno(:,:,end)); unique(anno(:,1,:)); ...
    unique(anno(:,end,:)); unique(anno(1,:,:)); unique(anno(end,:,:))]);

% Find size filter
rpAnno = regionprops(anno,'PixelIdxList','Area');
uannoBig = find([rpAnno.Area] > 2000);

% All synapse partners
synPartners = unique(eTest(:,1:2));

% Find synapses that are orphans or small
synOrphan = synPartners(~ismember(synPartners, uannoEdge) | ~ismember(synPartners, uannoBig));

synT = [];
synAll = (eTest(:,3));
for i = 1:length(synAll)
    i1 = find(eTest(i,1)==synOrphan);
    i2 = find(eTest(i,2)==synOrphan);
    if length(i1)+length(i2) > 0
        synT = [synT, synAll(i)];
    end
end


%% Find parents; find merges

% Need to do reverse lookup
pLookup = [];
for i = 1:length(rpAnno)
    if rpAnno(i).Area > 0
        pLookup(end+1,:) = [i, mode(annoTruth.data(rpAnno(i).PixelIdxList))];
    end
end

%% Do automerge of all orphans
annoMergeAll = annoDetect.data;
for i = 1:length(synOrphan)
    idx = find(synOrphan(i) == pLookup);
    parent = pLookup(idx,2);
    pids = find(parent == pLookup(:,2));
    
    mergeFrom = pLookup(idx,1);
    
    
    mergeTo = (pLookup(pids,1));
    %length(mergeFrom)
    clear sizeA
    for j = 1:length(mergeTo)
        if mergeFrom == mergeTo
            sizeA(j) = 0;
        else
            sizeA(j) = annoMergeAll(rpAnno(mergeTo(j)).Area);
        end
    end
    iidx = find(sizeA == max(sizeA));
    iidx = iidx(1);
    mergeTo = mergeTo(iidx);
    annoMergeAll(rpAnno(mergeFrom).PixelIdxList) = mergeTo;
end
ANNOMergeAll = RAMONVolume; ANNOMergeAll.setCutout(annoMergeAll);

eTest = statconn_associate_neusyn(ANNOMergeAll, synTruth);
[~,~, lgTestMergeAll, ~] = statconn_construct_graph(eTest);

% Compute Error
gErrMergeAll = statconn_compute_graph_error_simple(lgTestMergeAll, lgTruth);

%% Do merges automatically with Gala (filtered)

% Load Gala Graphs
f = dir('v29*.mat');
for i = 1:length(f)
    load(f(i).name);
    gg(:,:,i) = gErrMetrics.lgTest>0;
end

mgg = mean(gg,3);
figure, imagesc(mgg+mgg')

count = 1;
thresh = unique(mgg);
clear TPAll scoreF1
baseGraph = lgAnnoDetect;

for i = 1:length(thresh)
    gd = (mgg >= thresh(i));
    gtfb = single((lgTruth+lgTruth')>0);
    gdfb = single((gd+gd')>0);
    
    gdfbf = zeros(size(gdfb));
    gdfbf(:,synT) = gdfb(:,synT); %maybe assumes synIndexing
    gdfbf(synT,:) = gdfb(synT,:);
    
    gdfbf = gdfbf + baseGraph + baseGraph' > 0;
    TP = sum(sum(gdfbf == 1 & gtfb == 1));
    FP = sum(sum(gdfbf == 1 & gtfb == 0));
    FN = sum(sum(gdfbf == 0 & gtfb == 1));
    prec = TP/(TP + FP);
    rec = TP/(TP + FN);
    scoreF1(i) = 2 * (prec * rec) / (prec + rec);
    TPAll(i) = TP;
    P(i) = prec;
    R(i) = rec;
end

fprintf('F1 Score for best Gala operating point: %0.4f\n',gErrAnnoDetect.scoreF1)
fprintf('F1 Score for ALL spine/bouton merge operating point: %0.4f\n',gErrMergeAll.scoreF1)
fprintf('F1 Score for automatic spine/bouton merge operating point (based on Gala graphs): %0.4f\n',max(scoreF1))


figure(101)

subplot(2,2,1)
imagesc(lgTruth+lgTruth')
title('True Line Graph')

subplot(2,2,2)
imagesc(lgAnnoDetect+lgAnnoDetect')
title('Gala (0.55 Agglo Threshold) Line Graph [f1 score: 0.32]')

subplot(2,2,3)
imagesc(lgTestMergeAll+lgTestMergeAll')
title('Line Graph - Magically Merge Boutons and Spines [f1 score: 0.62]')

subplot(2,2,4)
imagesc(mgg+mgg')
title('Probabalistic (Average) Gala Line Graph')

set(gcf,'Color',[1 1 1])