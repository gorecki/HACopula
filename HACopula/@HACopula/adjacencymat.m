function adjacencyMatrix = adjacencymat(obj)
% ADJACENCYMAT - a d x d matrix showing 1 if two nodes (leaves or forks) are connected in
% HAC structure and 0, otherwise.
%
% Note that the forks are labeled according to TauOrdering.
%
% Example:
%
% myHAC = HACopula({{'A', .5}, 1, {{'A', 0.9}, 2, {{'C', 2}, 3, 4}}});
% adjacencymat(myHAC)
%
% ans =
%
%      0     0     0     0     0     0     1
%      0     0     0     0     0     1     0
%      0     0     0     0     1     0     0
%      0     0     0     0     1     0     0
%      0     0     1     1     0     1     0
%      0     1     0     0     1     0     1
%      1     0     0     0     0     1     0
%
%
% Copyright 2017 Jan Górecki and Martin Holeòa


orderedNodes = obj.Forks;
nNodes = size(orderedNodes, 2) + size(obj.Leaves,2);
adjacencyMatrix = zeros(nNodes, nNodes);
for node = 1:length(orderedNodes)
    thisNode = orderedNodes{node};
    childArray = zeros(1, size(thisNode.Child, 2));
    for j = 1:length(thisNode.Child)
        if isa(thisNode.Child{j}, 'HACopula')
            childArray(j) = thisNode.Child{j}.TauOrdering;
        else
            childArray(j) = thisNode.Child{j};
        end
    end
    
    adjacencyMatrix(thisNode.TauOrdering, childArray) = ones(1, size(childArray,2));
end
% make it symetric
adjacencyMatrix = adjacencyMatrix + adjacencyMatrix';
end