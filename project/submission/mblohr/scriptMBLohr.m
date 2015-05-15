%% scriptMBLohr.m
% My final project for statconn evaluates a connectivity brain graph 
% feature in epileptic patients, to determine if they will respond 
% positively (seizure or AD is mitigated) or negatively (seizure or AD
% is not mitigated) following a mental effort exercise. 
% Michele Lohr; May 14th, 2015
%
%The script must be run from in a directory structured as follows:
% **(w/ other files sure, no other directories)**
%
%-Current Directory:
%   -scriptMBLohr.m
%   [the following are contents of a ZIP this script will download]
%   -dataMBLohr.zip
%   -dataMBLohr.mat
% 
% ... only external dependency is wget

% %% Preamble
% 
% %have I downloaded my zip file
fname = 'dataMBLohr.zip';
files = dir();

clear found
for i=1:length(files) %silly but your directory should be almost empty
    if strcmp(fname, files(i).name)
        found = 1;
        break;
    end
end

%if no, download it
if ~exist('found', 'var')
    %craft and run command
    addr = ' https://www.dropbox.com/sh/idt3d0gylplyo31/AABrSIY6piow4kA64TUcAvFUa/dataMBLohr.zip';
    opts = ' --no-check-certificate';
    cmd = strcat('wget', addr, opts);
    system(cmd, '-echo');
end

%regardless, unzip it
unzip(fname);

%% Evaluate
load dataMBLohr  % loads electrodes, nonresponsive_EVCs, and responsive_EVCs

figure(1)
subplot(1,2,1)
imagesc(nonresponsive_EVCs');
xlabel('Event No.')
ylabel('Electrode No.') 
title('Non-Responsive')
subplot(1,2,2)
imagesc(responsive_EVCs');
xlabel('Event No.')
ylabel('Electrode No.') 
title('Responsive')
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,['\bf EVC Values'],'HorizontalAlignment','center','VerticalAlignment', 'top');
    drawnow

for j = 1:88
    for k = 1:88
        distMatrix(j,k) = sqrt((electrodes(j,1)-electrodes(k,1))^2 + (electrodes(j,2)-electrodes(k,2))^2);
    end
end

[numNonResponsive cols] = size(nonresponsive_EVCs);
[numResponsive cols] = size(responsive_EVCs);

for k = 1:numNonResponsive
    [Y,ind] = sort(nonresponsive_EVCs(k,:),'descend');
    for j = 1:19
        summation(j) = sum(distMatrix(ind(j),:));
    end
    totalDist_nonResponsive(k) = sum(summation);
end

for k = 1:numResponsive
    [Y,ind] = sort(responsive_EVCs(k,:),'descend');
    for j = 1:19
        summation(j) = sum(distMatrix(ind(j),:));
    end
    totalDist_Responsive(k) = sum(summation);
end

figure(2)
plot([1:numNonResponsive],totalDist_nonResponsive,'b.-',...
    [1:numResponsive],totalDist_Responsive,'r.-');hold on;
plot([0 12],[10050 10050],'k--');
text(6,10100,'Proposed Decision Boundary');
text(8,10000,'(Patient-Specific)')
xlabel('Event No.');
ylabel('Summed Distances Between High EVC Electrodes');
legend('Non-Responsive','Responsive');
title('\bf EVC Distance Feature');
   
