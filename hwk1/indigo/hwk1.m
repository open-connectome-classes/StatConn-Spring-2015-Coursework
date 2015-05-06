%The following is an example when kmeans "fails", meaning that it doesn't
%cluster this graph in a way that most humans would 

%If a human were clustering, they would probably place nodes A and B into a
%cluster, D, E and F into a cluster, and C could go in either cluster. This
%seems to most reasonable, yet kmeans fails to produce these results.
     

%    A B C D E F     
a = [0,1,1,0,0,0; %A
     1,0,1,0,0,0; %B
     1,1,0,1,1,1; %C
     0,0,1,0,1,1; %D
     0,0,0,1,0,1; %E
     0,0,0,1,1,0]; %F
 
%Run k means to get the clusters 
clusters = kmeans(a, 2)

%The output always places nodes A, B, and D in a cluster and C, E, and F in
%a cluster. This seems strange, and it would be counterintuative for a
%human to cluster in this way.