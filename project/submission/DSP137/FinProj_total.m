% (1) Kernel Based Hypothesis Testing: Maximum Mean Discrepancy
% Using T_n = (n1+n2)(norm(mu2-mu1)_H)^2 as test statistic - MMD
% (1i) Have sum of BOLD levels for each region in each patient:
% Proceed to normalize BOLD levels across patient. Use this as a way to get
%   an 'average' region i measure for person j
% (1ii) Have average of BOLD values across region for each participant.
%
% HGP 20150508
% S.D.G.

url='https://www.dropbox.com/sh/idt3d0gylplyo31/AABON28Y7zVAqWK_dk3KN7XSa/DSP137_heathergaddypatsolic_finproj694.zip?dl=1';
cmd = ['wget ' url ' -O "DSP137_heathergaddypatsolic_finproj694.zip" --no-check-certificate'];
system(cmd, '-echo');
unzip('DSP137_heathergaddypatsolic_finproj694.zip');
addpath('DSP137_heathergaddypatsolic_finproj694');

load('partroisum_matrix.mat'); % This loads matrix vecs
load('avgvecs.mat'); % This loads avgvecs matrix
%% Create Label Vector for ADHD v. TDC
M = csvread('ADHDlabels_ordered.csv',1,0);
M = M(:,[1,3]);

v=M(:,1);
I = eye(length(v));
P = I(v,:);
M2 = P'*M;
PRmatrixO = P'*vecs; %particpant roi matrix ordered 1 to 40


% Re-Ordered matrix so that ADHD on top and TDC on bottom +
% ADHD row i matched with TDC row 20+i
newmat = zeros(size(vecs));
for i = 1:40
    if M2(i,2) == 1
        j = (i+1)/2;
        newmat(j,:) = PRmatrixO(i,:);
    else %M2(i,2) == 0
        j = 20 + (i/2);
        newmat(j,:) = PRmatrixO(i,:);
    end
end

%% Do same as above for the avgvecs matrix
PRmatrixOavg = P'*avgvecs; %particpant roi matrix ordered 1 to 40


% Re-Ordered matrix so that ADHD on top and TDC on bottom +
% ADHD row i matched with TDC row 20+i
newmatavg = zeros(size(avgvecs));
for i = 1:40
    if M2(i,2) == 1
        j = (i+1)/2;
        newmatavg(j,:) = PRmatrixOavg(i,:);
    else %M2(i,2) == 0
        j = 20 + (i/2);
        newmatavg(j,:) = PRmatrixOavg(i,:);
    end
end

%%

% normalized matrix for newmat. 
% normalize the row corresponding to each patient.
nnewmat = zeros(size(newmat));
for i=1:40
    normi = norm(newmat(i,:));
    nnewmat(i,:) = newmat(i,:)/normi;
end

%% Calculate MMD on nnewmat and run hypothesis testing 
%  p-val plot should match first figure in the report
%  For each region:
T = zeros(70,1); % vector of gaussian kernel estimates for each region
%no = zeros(70,1); % will be a vector of difference norms
pval = zeros(70,1);
for r = 1:70
    if nnewmat(:,r) == zeros(40,1)
        pval(r) = 1;
    else
        adhd = nnewmat(1:20,r);
        tdc = nnewmat(21:40,r);
        gamma = 1/100; % this can be changed, note: gamma=1/(2(\sigma^2))
        n = norm(adhd-tdc);
        %no(r) = n;
        %% Compute similarity factor using Gaussian kernel
        gkern = exp(-gamma*n^2); % This is the Hilbert Space Kernel being used
        T(r) = 40*gkern^2; % This is the test statistic
        %% permutation test
        Tperm = [];
        for ii = 1:10000
            N = randperm(40);        
            P2 = I(N,:);
            PRmatrixOperm = P2'*nnewmat;
            adhdperm = PRmatrixOperm(1:20,r);
            tdcperm = PRmatrixOperm(21:40,r);
            gk = exp(-gamma*norm(adhdperm-tdcperm)^2);
            t = 40*gk^2;
            Tperm = [Tperm; t];         
        end
        
        pval(r) = sum(Tperm >= T(r)) ./ length(Tperm);
        
        %% If desired, can print out kernel density plot for the permutation
        %     T statistics to see the curve and test statistic for each ROI
        %         ksdensity(Tperm)
        %         hold on
        %         plot(T(r),0,'x')
        %         title(strcat('Kernel Density plot for Test Statistic for ROI: ',num2str(r),'(normalized version)'));
        %         hold off
        %         pause
    end
end

figure%(2)
hold on
plot(pval)
for r=1:70
    if pval(r) < .01
        plot(r,pval(r),'x')
    end
end
title(strcat('MMD test statistic P-value plot'));
xlabel('Region of Interest');
ylabel('p-value');
hold off


%% Calculate MMD on newmatavg and run hypothesis testing 
%  p-val plot should match second figure in the report
Tavg = zeros(70,1); % vector of gaussian kernel estimates for each region
%no = zeros(70,1); % will be a vector of difference norms
pval = zeros(70,1);
for r = 1:70
    if newmatavg(:,r) == zeros(40,1)
        pval(r) = 1;
    else
        adhd = newmatavg(1:20,r);
        tdc = newmatavg(21:40,r);
        gamma = 1/100; % this can be changed
        n = norm(adhd-tdc);
        %no(r) = n;
        %% Compute similarity factor using Gaussian kernel
        gkern = exp(-gamma*n^2); % This is the Hilbert Space Kernel being used
        Tavg(r) = 40*gkern^2; % This is the test statistic
        %% permutation test
        Tperm = [];
        for ii = 1:10000
            N = randperm(40);
            %tempg = gErr(N, N); <- not sure what this line does
            %     figure(4); imagesc(tempg)
            
            P2 = I(N,:);
            PRmatrixOperm = P2'*newmatavg;
            adhdperm = PRmatrixOperm(1:20,r);
            tdcperm = PRmatrixOperm(21:40,r);
            gk = exp(-gamma*norm(adhdperm-tdcperm)^2);
            t = 40*gk^2;
            Tperm = [Tperm; t];
            set(gca,'XLimMode','manual','YLimMode','manual');  % Fix axes limits
            
        end
        
        pval(r) = sum(Tperm >= Tavg(r)) ./ length(Tperm);
        
%         ksdensity(Tperm)
%         hold on
%         plot(T(r),0,'x')
%         title(strcat('Kernel Density plot for Test Statistic for ROI: ',num2str(r),'(normalized version)'));
%         hold off
%         pause
    end
end 

figure%(3)
hold on
plot(pval)
for r=1:70
    if pval(r) < .01
        plot(r,pval(r),'x')
    end
end
title(strcat('MMD test statistic for Average BOLD across ROI: P-value plot'));
xlabel('Region of Interest');
ylabel('p-value');
hold off

%% Hellinger Distance - The code for this section was modified from {
% W Gray Roncal
% 01.14.2015
% Greg Kiar
% 05.04.2015
% Script to compare KKI42 test retest data and plot result
% Scans are reordered so that both scans for each subject appear together
% Data (covariates and small graphs are from the public ftp site)
% openconnecto.me/data/public/MR/MIGRAINE_v1_0/KKI-42/
% }
% https://github.com/openconnectome/m2g/blob/master/test/analysis/kki_compare.m

pval = zeros(70,1);
for r = 1:70; 
    if newmat(:,r) == 0;
        pval(r) = 1;
    else
        adhd = newmat(1:20,r);
        tdc = newmat(21:40,r);
        
        %figure(1)
        [~, x_adhd] = ksdensity(adhd);
        [~, x_tdc] = ksdensity(tdc);
        lims = [min([x_adhd,x_tdc]), max([x_adhd, x_tdc])]; %innefficient but works...
        xrange = lims(1):range(lims)/300:lims(2);
        
        [f_adhd] = ksdensity(adhd, xrange);
        [f_tdc] = ksdensity(tdc, xrange);
        
        H = norm(sqrt(f_adhd)- sqrt(f_tdc),2)/sqrt(2); %computed Hellinger distance
        
%         % For visual if desire to see comparison of kdensity for TDC v
%         % ADHD per ROI
%         plot(xrange, f_adhd, xrange, f_tdc);
%         title(strcat('Hellinger Distance=', num2str(H),'  for ROI: ',num2str(r)));
%         legend('ADHD Subject Kernel Esimate', 'TDC Subject Kernel Estimate');
%         pause
        
        %% Compute Null Distribution and Hellinger
        
        hperm = [];
        for ii = 1:1000
            N = randperm(40);
            P1 = I(N,:);
            newmatperm = P1'*newmat; %permute the rows (still use top 20 as ADHD)
            
            n_adhd = newmatperm(1:20,r);
            n_tdc = newmatperm(21:40,r);
            
            %        figure(2)
            [~, x1_adhd] = ksdensity(n_adhd);
            [~, x1_tdc] = ksdensity(n_tdc);
            lims = [min([x1_adhd,x1_tdc]), max([x1_adhd, x1_tdc])]; %innefficient but works...
            nxrange = lims(1):range(lims)/300:lims(2);
            
            [f_n_adhd] = ksdensity(n_adhd, nxrange);
            [f_n_tdc] = ksdensity(n_tdc, nxrange);
            
            % Compute Hellinger distance for each permuted labels matrix and
            % add to hperm storage vector
            hperm = [hperm, norm(sqrt(f_n_adhd)- sqrt(f_n_tdc),2)/sqrt(2)];
            %         if rem(ii,100)==0
            %             figure(2); plot(nxrange, f_n_adhd, nxrange, f_n_tdc);
            %             title(strcat('Perm: Hellinger Distance=', num2str(H),'  for ROI: ',num2str(r)));
            %             legend('ADHD Subject Kernel Esimate', 'TDC Subject Kernel Estimate');
            %             pause
            %         end
        end
        pval(r) = sum(hperm >= H) ./ length(hperm);
    end
end
figure%(4)
hold on
plot(pval)
for r=1:70
    if pval(r) < .01
        plot(r,pval(r),'x')
    end
end
title(strcat('Hellinger Distance P-value plot for each ROI'));
xlabel('Region of Interest');
ylabel('p-value');
hold off

%% Wilcoxon Signed-Rank Test for paired samples
% For each region:
% compute absolute value vector for each region matrix is 20x70
% each column represents the vector of matched differences for that region
% similar for sgn matrix (20x70)
% absvalm = zeros(20,70);
% sgm = zeros(20,70);
Wsgnrankpval = zeros(70,1); %holds pvalue for each region of interest
for r = 1:70
    adhd = newmat(1:20,r);
    tdc = newmat(21:40,r);
    Wsgnrankpval(r) = signrank(adhd,tdc);
end

figure%(5)
hold on
plot(Wsgnrankpval)
for r=1:70
    if Wsgnrankpval(r) < .01
        plot(r,pval(r),'x')
    elseif Wsgnrankpval(r) < 0.05 && Wsgnrankpval(r)>0.01
        plot(r,Wsgnrankpval(r),'o')
    end
end
title(strcat('Wilcoxon Signed-Rank Test P-value plot for each ROI'));
xlabel('Region of Interest');
ylabel('p-value');
hold off



