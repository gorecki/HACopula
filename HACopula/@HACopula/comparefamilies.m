function [areEqual, ratio] = comparefamilies(HAC1, HAC2)
%COMPAREFAMILIES - compare the family of corresponding forks in two HACs
%
% If the structures of two HACopula inputs HAC1 and HAC2 are the same
% (i.e., comparestructures(HAC1, HAC2) = 1), then the method for each fork
% in HAC1 compares the family of this fork to the family of the
% corresponding fork (according to the structure) in HAC2. Then, it returns
% areEqual = 1, if all of them are the same, areEqual = 0, otherwise.
% If the second output (ratio) is required, it returns the ratio of forks
% with the same families to all forks.
%
%
% Copyright 2018 Jan Gorecki

if ~comparestructures(HAC1, HAC2)
    error('HACopula:BadInputs', 'HACopula::comparefamilies: HAC1 and HAC2 have to be of the same structure (can be checked by the method ''comparestructures'').');
end
% now HAC1.Dim = HAC2.Dim and also both HACs have the same number of forks

d = HAC1.Dim;
k = size(HAC1.Forks,2);
leaves1 = cell(1, k);
leaves1Extend = zeros(k, d);
families1 = cell(1, k);
leaves2 = cell(1, k);
leaves2Extend = zeros(k, d);
families2 = cell(1, k);
for i = 1:k
    leaves1{i} = HAC1.Forks{i}.Leaves;
    leaves2{i} = HAC2.Forks{i}.Leaves;
    
    % extend leaves for sortrow below
    leaves1Extend(i, 1:d) =  [zeros(1, d - size(leaves1{i},2)) sort(leaves1{i})];
    leaves2Extend(i, 1:d) =  [zeros(1, d - size(leaves2{i},2)) sort(leaves2{i})];
    
    % get the families
    families1{i} = HAC1.Forks{i}.Family;
    families2{i} = HAC2.Forks{i}.Family;
end

[~, index1] = sortrows(leaves1Extend);
[~, index2] = sortrows(leaves2Extend);

% sort them
sortedFamilies1 = families1(index1);
sortedFamilies2 = families2(index2);
ratio = sum(strcmp(sortedFamilies1, sortedFamilies2))/k;
areEqual = (ratio == 1);

end
