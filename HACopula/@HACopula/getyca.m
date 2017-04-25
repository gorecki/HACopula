function [ycaTauOrdering, ycaFork] = getyca(obj, leaf1, leaf2)
%GETYCA - returns in ycaTauOrdering the TauOrdering property of the youngest common
% ancestor of the leaves leaf1 and leaf2.
% Also the youngest common ancestor fork is returned in ycaFork.
%
%
% Copyright 2017 Jan Górecki

if (size(intersect([leaf1 leaf2], obj.Leaves), 2) ~= 2)
    error('HACopula::getycaindex: leaf1 and leaf2 must be two different leaves of the HAC.');
end
orderedHACModel = obj.Forks;
d = getdimension(obj);
k = size(orderedHACModel, 2);
minLeavesSize = d+1;
% get all sets of the leaves
for i = 1:k
    leavesVec = orderedHACModel{i}.Leaves;
    if (size(intersect([leaf1 leaf2], leavesVec), 2) == 2) && ...
            (size(leavesVec, 2) < minLeavesSize)
        % yca is the fork with the fewest descendat leaves
        % containing both leaf1 and leaf2
        ycaTauOrdering = orderedHACModel{i}.TauOrdering;
        ycaFork = orderedHACModel{i};
    end
end
end