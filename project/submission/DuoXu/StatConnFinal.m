%% Import Data

clear
close all

%have I downloaded my zip file
fname = 'dx_statconn.zip';
files = dir();

for i=1:length(files) %silly but your directory should be almost empty
    if strcmp(fname, files(i).name)
        found = 1;
        break;
    end
end

%if no, download it
if ~exist('found', 'var')
    %craft and run command
    addr = ' https://www.dropbox.com/sh/idt3d0gylplyo31/AADwB7_QznoAKw9wFagFRW9Pa/dx_statconn.zip?dl=1';
    opts = ' --no-check-certificate';
    cmd = strcat('wget', addr, opts);
    system(cmd, '-echo');
end

%regardless, unzip it
unzip(fname);



graphNum = length(pp.btypeInd);
bTypeNames = { 'Hit', 'Miss', 'False Alarm', 'Correct Rejection'};
bTypeAbr = { 'H', 'M', 'FA', 'CR'};

set(0, 'DefaultFigureWindowStyle', 'docked');

figure(1)
pp.uSession.ShowCorrBin();
title('Functional (Pearson''s correlation) connectivity of neurons for all trials');

pause(.5);

figure(2)
for i = 1 : graphNum
    subplot(2,2,i);
    pp.uBtypes{i}.ShowCorrBin();
    title([ 'Functional connectivity in ' bTypeNames{i} ' trials' ]);
end

pause(.5);


%% Fit SBM with Different Ks to Data

sbmKs = num2cell((1:8)');
sbmLabels = cell(length(sbmKs),1);
for i = 1 : length(sbmKs)
    sbmLabels{i} = kmeans(pp.uSession.pCorrBin, sbmKs{i});
end

sbms = cell(length(sbmKs), graphNum);
for k = 1 : length(sbmKs)   % for every cluster number, k
    for gDraw = 1 : graphNum   % for every behavioral types
        sbms{k,gDraw} = SBM(pp.uBtypes{gDraw}.pCorrBin);
        sbms{k,gDraw}.SetLabels(sbmLabels{k});
        sbms{k,gDraw}.FitData();
    end
end

figure(3)
for i = 1 : length(sbms(:))
    subplot(size(sbms,2), size(sbms,1), i);
    sbms{i}.ShowBernBlock();
    
    if ~mod(i-1, length(sbmKs))
        ylabel(bTypeNames{(i-1)/length(sbmKs)+1});
    end
    title([ 'k = ' num2str(sbms{i}.k) ]);
    set(gca, 'xTickLabel', {});
    set(gca, 'yTickLabel', {});
end

pause(.5);


%% Find Optimal K

logLik = cell(length(sbmKs),1);
for k = 1 : length(sbmKs)
    logLik{k} = zeros(graphNum);
    for i = 1 : graphNum
        for j = 1 : graphNum
            if i ~= j
                [ logLik{k}(i,j), aic{k}(i,j) ] = LogLik(sbms{k,j}, sbms{k,i}.adjMatrix);
            end
        end
    end
end

sumLogLik = cellfun(@(x) sum(x(:)), logLik);
[ sbmOptKVal, sbmOptKInd ] = max(sumLogLik);
sbmOptK = sbmKs{sbmOptKInd};

figure(4)
plot(sumLogLik);
hold on
plot(sbmOptKInd, sbmOptKVal, 'o');
text(sbmOptKInd+0.2, sbmOptKVal, [ 'Optimal k = ' num2str(sbmOptKInd) ]);
hold off

xlabel('k');
ylabel('\Sigma\Sigma logP (G(i) | SBM(j,k)), i ~= j');
title('Sum Log-likelihood of Cross-sampling w.r. Cluster Number k');

pause(.1);


%% Parametric Bootstrap with K*

%{
Fij = log(f(Gi|theta(k*)(Gj)))
Pij = p-value
%}

bsN = 1000;
bsLogLiks = zeros(graphNum, graphNum, bsN);
bsPs = zeros(graphNum);

for gDraw = 1 : graphNum
    for gTest = 1 : graphNum
        for s = 1 : bsN
            bsSample = sbms{sbmOptKInd,gDraw}.DrawAdjMat();
            bsLogLiks(gDraw,gTest,s) = LogLik(sbms{sbmOptKInd,gTest}, bsSample);
        end
        bsLogLiks(gDraw,gTest,:) = sort(bsLogLiks(gDraw,gTest,:));
    end
end

figure(5)
% Given 'Title' modeled with k*, the log-probability of 
for gDraw = 1 : graphNum
    subplot(4, 3, [ gDraw*3-2 gDraw*3-1 ]);
    for gTest = 1 : graphNum
    hold on
        h = histogram(squeeze(bsLogLiks(gDraw,gTest,:)));
        text(mean(h.BinLimits), max(h.Values)*1.05, [ '#' num2str(gTest) ]);
    end
    hold off
    xlabel([ 'logP (G(#) | SBM(' bTypeAbr{gDraw} ',k=' num2str(sbmOptKInd) ')' ]);
    ylabel('Count');
    
    subplot(4, 3, gDraw*3);
    text(0, 0.5, [ '#' num2str(gDraw) ': ' bTypeNames{gDraw} ' (' bTypeAbr{gDraw} ')' ], 'FontSize', 12);
    axis off
end


