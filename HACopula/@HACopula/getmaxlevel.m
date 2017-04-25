function maxLevel = getmaxlevel(obj)
%GETMAXLEVEL - the maximal nesting level in obj.
%
%
% Copyright 2017 Jan Górecki

maxLevel = -Inf;
for i = 1:size(obj.Forks,2)
    maxLevel = max(maxLevel, obj.Forks{i}.Level);
end
end