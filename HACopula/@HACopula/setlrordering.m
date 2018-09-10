function setlrordering(obj)
%SETLRORDERING - switch children of a fork in order to canonize
%                the graphical representation of a HAC.
%
% This auxiliary function (for the plot function) interchanges
% the children of each fork in order to get a left-right
% ordering, i.e., a child with less dencendant forks are placed
% at the left-side of a child with more dencendant forks. If
% the numbers of the dencendant forks are equal for two or more
% children, the child with the lowest descendant leaf index is
% placed at the left-side of the remaining children (which are
% ordered analogously).
%
% This ordering assures that a fully-nested Archimedean copula
% is plotted (by the plot function) as is common on "the main
% diagonal"
%
% NOTE: This function modifies the input obj.
%
%
% Copyright 2018 Jan Gorecki

nChildren = size(obj.Child, 2);
nDescForks = zeros(1, nChildren);
lowestLeafIndex = zeros(1, nChildren);
for i = 1:nChildren
    if isa(obj.Child{i},'HACopula')
        nDescForks(i) = size(obj.Child{i}.Forks,2) - 1;
        lowestLeafIndex(i) = min(obj.Child{i}.Leaves);
    else
        nDescForks(i) = -1;
        lowestLeafIndex(i) = obj.Child{i};
    end
end
[~, lrOrdering] = sortrows([nDescForks' lowestLeafIndex']);
% interchange the children according to lrOrdering
child = cell(1, nChildren);
for i = 1:nChildren
    child{i} = obj.Child{lrOrdering(i)};
end
obj.Child = child;

% do recursion
for i = 1:nChildren
    if isa(obj.Child{i},'HACopula')
        setlrordering(obj.Child{i});
    end
end
end