function [ycaTauOrdering, ycaFork] = getyca(obj, leaf1, leaf2)
%GETYCA - returns in ycaTauOrdering the TauOrdering property of the youngest common
% ancestor of the leaves leaf1 and leaf2.
% Also the youngest common ancestor fork is returned in ycaFork.
%
% Copyright 2018 Jan Gorecki

if not(any(leaf1 == obj.Leaves) && any(leaf2 == obj.Leaves) && leaf1 ~= leaf2)
    error('HACopula:getycaindex', 'HACopula::getycaindex: leaf1 and leaf2 must be two different leaves of the HAC.');
end

orderedHACModel = obj.Forks;
d = obj.Dim;
k = size(orderedHACModel, 2);
minLeavesSize = d+1;
% get all sets of the leaves
for i = 1:k
    leavesVec = orderedHACModel{i}.Leaves;
    if any(leaf1 == leavesVec) && any(leaf2 == leavesVec) && ...
       (size(leavesVec, 2) < minLeavesSize)
        % yca is the fork with the fewest descendat leaves
        % containing both leaf1 and leaf2
        ycaTauOrdering = orderedHACModel{i}.TauOrdering;
        ycaFork = orderedHACModel{i};
        minLeavesSize = size(leavesVec, 2);
    end
end

end