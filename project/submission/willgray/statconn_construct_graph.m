function [neuGraph, nId, synGraph, synId] = statconn_construct_graph(edgeList, varargin)

if size(edgeList,2) == 5 && sum(edgeList(:,5)) ~= 0
    eField = 5;
else
    eField = 3;
end
synId = unique(edgeList(:,eField));

synId(synId == 0) = [];

synGraph = zeros(length(synId));

neuConn = edgeList(:,1:2);
nId = unique(neuConn);
for i = 1:length(nId)
    
    [dataR, ~] = find(neuConn == nId(i));
    
    sIdMatch = unique(edgeList(dataR,eField));
    
    combo = combnk(sIdMatch,2);
    combo = sort(combo,2);
    for j = 1:size(combo,1)
        s1 = find(synId == combo(j,1));
        s2 = find(synId == combo(j,2));
        synGraph(s1,s2) = synGraph(s1,s2)+1;
    end
end

%All neurons not represented have degree 0
% syn-neuron: Any neuron
neuGraph = zeros(length(nId));
for i = 1:size(neuConn,1)
    n1 = neuConn(i,1);
    n2 = neuConn(i,2);
    
    n1 = find(n1 == nId);
    n2 = find(n2 == nId);
    neuGraph(min(n1,n2),max(n1,n2)) = neuGraph(min(n1,n2),max(n1,n2)) + 1;
end
