function decreaselevel(obj)
%DECREASELEVEL - decreases the level of all forks in obj.
%
% Used after collapsing of two forks in the method collapse.
%
% NOTE: This method modifies the input obj.
%
%
% Copyright 2017 Jan Górecki

for i = 1:size(obj.Forks,2)
    if isa(obj.Forks{i},'HACopula')
        obj.Forks{i}.Level = obj.Forks{i}.Level - 1;
    end
end
end