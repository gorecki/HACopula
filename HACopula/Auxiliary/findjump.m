function jumpIndex = findjump(distArray)
%FINDJUMP - Find the first substantial jump in the minimal distances
%
% Purpose:
% This function is used to heuristically find an appropriate number of
% forks while collapsing a HAC. More precisely, based on distArray
% (obtained using [~, distArray] = collapse(obj ,...)), jumpIndex =
% findjump(distArray) returns the index of the first step
% (=distArray(k+1)-distArray(k)) that is higher than the average step
% (=distArray(1)-distArray(end)/length(distArray)). For details, see
% Section 6.1 in 
%
% References:
% [Gorecki et al., 2017] On Structure, Family and Parameter Estimation
%     of Hierarchical Archimedean copulas. Journal of Statistical Computation 
%     and Simulation, 87(17), 3261-3324
%
%
% Copyright 2018 Jan Gorecki

if size(distArray,2) == 1
    jumpIndex = 1;
    return
end

firstInd = 1;
if distArray(1) == 0
    % find the first value higher than 0 
    firstInd = find(distArray > 0 , 1, 'first');
    % and reduce the array
    distArrayOrig = distArray;
    distArray = distArray(firstInd:end);
end
if size(distArray,2) > 0
    first = distArray(1);
else % there are no positives in distArray
    jumpIndex = size(distArrayOrig,2); % collapse to maximum
    return
end

n = length(distArray);
last = distArray(n);

    
MIN_OVERALL_RATIO = 0.5;
if (last - first)/first < MIN_OVERALL_RATIO % if the difference is very small
    % it looks like that it is an AC, i.e., there is only one fork
    jumpIndex = n;
else
    step = (last - first)/(n-1);
    diffDistArray = distArray(2:n)-distArray(1:n-1);
    jumpIndex = find(diffDistArray >= step , 1, 'first');
end

% transform to the original size
jumpIndex = jumpIndex + firstInd - 1;

end