function [areEqual, areEqualFamilies, nEqualFamilies] = comparestructures(HACObject1, HACObject2)
%COMPARESTRUCTURES - returns areEqual = 1, if the structures of HACObject1 and
% HACObject2 are equal, otherwise, returns areEqual = 0.
%
% The comparison is based on comparing the sets of the
% descendant leaves of all forks. Additionally, if the
% structures are equal, the method compares the families of the
% corresponding generators, and returns areEqualFamilies = 1,
% if the families are all equal, areEqualFamilies = 0, at least
% one family is different, nEqualFamilies is the number of the
% equal families. If the structures are not equal,
% areEqualFamilies = nEqualFamilies = NaN.
%
%
% Copyright 2017 Jan Górecki

d1 = getdimension(HACObject1);
d2 = getdimension(HACObject2);

areEqualFamilies = NaN;
nEqualFamilies = NaN;

if (d1 ~= d2)
    areEqual = false;
    return;
end

d = d1;

forks1 = HACObject1.Forks;
forks2 = HACObject2.Forks;

if (size(forks1,2) ~= size(forks2,2))
    areEqual = false;
    return;
end

k = size(forks1,2);

% get leaves of each of the forks1
leaves1 = cell(1, k);
leaves1Extend = zeros(k, d);
families1 = cell(1, k);
leaves2 = cell(1, k);
leaves2Extend = zeros(k, d);
families2 = cell(1, k);
for i = 1:size(forks1,2)
    leaves1{i} = forks1{i}.Leaves;
    leaves2{i} = forks2{i}.Leaves;
    
    % extend leaves for sortrow below
    leaves1Extend(i, 1:d) =  [zeros(1, d - size(leaves1{i},2)) sort(leaves1{i})];
    leaves2Extend(i, 1:d) =  [zeros(1, d - size(leaves2{i},2)) sort(leaves2{i})];
    
    % get the families
    families1{i} = forks1{i}.Family;
    families2{i} = forks2{i}.Family;
end

[leaves1Sorted, index1] = sortrows(leaves1Extend);
[leaves2Sorted, index2] = sortrows(leaves2Extend);
areEqual = (sum(sum(leaves1Sorted == leaves2Sorted)) == k*d);
if areEqual
    sortedFamilies1 = families1(index1);
    sortedFamilies2 = families2(index2);
    nEqualFamilies = sum(strcmp(sortedFamilies1, sortedFamilies2));
    areEqualFamilies = (nEqualFamilies == k);
else
    nEqualFamilies = NaN;
    areEqualFamilies = NaN;
end
end