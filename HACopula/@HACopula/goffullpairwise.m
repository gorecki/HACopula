function Sn = goffullpairwise(obj, U, aggFcn, statisticName, emp2copulas)
%GOFFULLPAIRWISE - Aggregated pairwise goodness-of-fit
%statistic
%
% This function computes the bivariate version of the Cramér-von Mises
% statistic statisticName ('K', 'R' or 'E') for each of the
% bivariate margins of the HAC obj and the observations U, and
% aggregates them using the aggregation function aggFcn
% ('average', 'max' or 'mean').
%
% NOTE:
% The input emp2copulas is optional and can be pre-computed using
% computeallemp2copulas(U).
%
%
% Copyright 2017 Jan Górecki

% perform basic data checks
if ~iscopuladata(U)
    warning('HACopulafit::goffullpairwise: At least one column of U was rejected to be standard uniform according to KS test at the 5% significance level, i.e., U might not be a sample from a copula.');
end

if (~exist('emp2copulas', 'var'))
    emp2copulas = computeallemp2copulas(U);
end

d = size(U, 2);
if getdimension(obj) ~= d
    error('HACopula::goffullpairwise: The data and the HAC are of different dimensions.')
end

pairSn = zeros(d, d);
for i = 1:d
    for j = i+1:d
        fork = getforkwithtauindex(obj, getyca(obj, i, j));
        pairSn(i, j) = gof2Sn(U(:, [i j]), fork.Family, fork.Parameter, statisticName, emp2copulas{i,j});
    end
end

% get upper triangular part of pairSn (without the main diagonal)
B = pairSn';
pairSnVec = B(tril(true(size(B)),-1));

switch aggFcn
    case 'max'
        Sn = max(pairSnVec);
    case 'average'
        Sn = mean(pairSnVec);
    case 'min'
        Sn = min(pairSnVec);
    otherwise
        error('gof2Sn: Unsupported aggregation function.');
end
end