function fork = getforkwithtauindex(obj, tauIndex)
%GETFORKWITHTAUINDEX - returns the fork with TauOrdering equal to tauIndex.
%
%
% Copyright 2017 Jan Górecki

for i = 1:length(obj.Forks)
    if obj.Forks{i}.TauOrdering == tauIndex
        fork = obj.Forks{i};
        return
    end
end
end