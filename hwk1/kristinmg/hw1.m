%% Statistical Connectomics
% Homework 1
% Kristin Gunnarsdottir

% Creating adjacency matrix with 3 clusters
A = zeros(10);
A(1,2) = 1;
A(2,1) = 1;
A(1,3) = 1;
A(3,1) = 1;
A(2,3) = 1;
A(3,2) = 1;
A(4,5 ) = 1;
A(5,4) = 1;
A(5,6) = 1;
A(6,5) = 1;
A(7,5) = 1;
A(5,7) = 1;
A(4,7) = 1;
A(7,4) = 1;
A(8,9) = 1;
A(9,8) = 1;
A(8,10) = 1;
A(10,8) = 1;

% A =
% 
%      0     1     1     0     0     0     0     0     0     0
%      1     0     1     0     0     0     0     0     0     0
%      1     1     0     0     0     0     0     0     0     0
%      0     0     0     0     1     0     1     0     0     0
%      0     0     0     1     0     1     1     0     0     0
%      0     0     0     0     1     0     0     0     0     0
%      0     0     0     1     1     0     0     0     0     0
%      0     0     0     0     0     0     0     0     1     1
%      0     0     0     0     0     0     0     1     0     0
%      0     0     0     0     0     0     0     1     0     0


% Using kmeans to partition the matrix into 3 clusters
kmeans(A,3)
% Here it gives me correct results:
% ans =
% 
%      1
%      1
%      1
%      2
%      2
%      2
%      2
%      3
%      3
%      3



kmeans(A,3)
% Here it fails:
% ans =
% 
%      1
%      1
%      1
%      2
%      2
%      2
%      2
%      3
%      1
%      1
%      

kmeans(A,3)
% Correct results:
% ans =
% 
%      2
%      2
%      2
%      1
%      1
%      1
%      1
%      3
%      3
%      3

kmeans(A,3)
% Fail:
% ans =
% 
%      1
%      1
%      1
%      2
%      3
%      2
%      2
%      1
%      1
%      1

