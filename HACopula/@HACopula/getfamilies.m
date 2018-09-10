function families = getfamilies(obj)
%GETFAMILIES - the set of families involved in the obj.
%
%
% Copyright 2018 Jan Gorecki

families = {};
for i = 1:length(obj.Forks)
    family = obj.Forks{i}.Family;
    if sum(strcmp(families, family)) == 0 % the family is not yet in families
        % add this family
        families = [{family} families];
    end
end
end