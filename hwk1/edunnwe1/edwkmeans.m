%% Kmeans assignment: Erika Dunn-Weiss
%
% Rather than generate a random graph, I drew a graph on paper that I
% expected to fail kmeans. Since this graph has two separated parts, I say
% that a node should be in cluster s if and only if it shares an
% edge with some other node within cluster s.
% I expected this graph to fail because the two clusters have large
% variance from the center of the cluster, which is an issue for kmeans
% since it is based on minimizing euclidean distance. 


load A
load Coordinates
% the graph in question:
figure; gplot(A, Coordinates,'-*')
ylim([0 6]), xlim([0 6]);
%%
IDX = kmeans(A, 2);
%%
% Given the graph, I expect the two disconnected parts to form each
% cluster. 
expected = vertcat(ones(6,1),2*ones(6,1)); %first version of expected results
expected2 = vertcat(2*ones(6,1),ones(6,1)); %second version of expected results
% test to see if clustered correctly
strcmp(IDX, expected)
strcmp(IDX, expected2)

%for me, comes out to 0 (false) every time. In fact, it comes out
%completely differently each time. Though, I could only run it a dozen
%times: for some reason when I try to put it in a for or while loop for a
%hundred iterations (on my matlab version 2012 I think) then it doesn't
% %work, with the error "??? Error using ==> kmeans>batchUpdate at 435
% Empty cluster created at iteration 1.
% 
% Error in ==> kmeans at 336
%         converged = batchUpdate();
% 
% Error in ==> edwkmeans at 19
% IDX = kmeans(A, 2);
% not sure what that's about. Anyway hope this suffices!
%%
% This section is new since our class discussion. I suspected that kmeans
% would still fail even if I initialized the clusters manually. This was
% indeed the case. I chose the two peaks of the triangle and the two
% centers of the line for two trials of manual initialization.
IDX1 = kmeans(A,2,'start',A([6,7],:));
strcmp(IDX1, expected)
strcmp(IDX1, expected2)

IDX2 = kmeans(A,2,'start',A([3,10],:));
strcmp(IDX2, expected)
strcmp(IDX2, expected2)