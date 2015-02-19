
numClasses = 2;

A_works = [0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 0 0 0; ...
           0 0 0 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
A_fails = [0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 0 0 1; ...
          0 0 0 0 1 1; 0 0 0 1 0 1; 0 0 1 1 1 0];
A = A_fails;

xy = [0 2; 1 2; 0.5 1; ...
      0 -2; 1 -2; 0.5 -1];

colorStr = {'b','r','g','k','y'};
colorStr2 = {'b.','r.','g.','k.','y.'};
idx = kmeans(A,numClasses);

figure(1);
subplot(1,numClasses+1,1);
gplot(A,xy)
hold on; plot(xy(:,1),xy(:,2),'b.');
title('Original Adjacency Matrix');
axis([-1 2 -3 3])
axis square;
axis off;

for ii = 1:numClasses
    subplot(1,numClasses+1,ii+1);
    idxClass = find(idx == ii);
    A2 = zeros(size(A));
    A2(idxClass,:) = A(idxClass,:);
    plot(xy(idxClass,1),xy(idxClass,2),colorStr2{ii});
    hold on;
    gplot(A2,xy,colorStr{ii});
    title(['Class ' num2str(ii)]);
    axis([-1 2 -3 3])
    axis square;
    axis off;
end



