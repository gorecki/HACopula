function out = findvars2flip(K)
%FINDVARS2FLIP - Find variables to flip in order to reduce negative pairwise correlations.
%
% Returs an array of indices of variables, which should be flipped (U := 1
% -U) in order to reduce negative correlation in data corresponding to the
% Kendall correlation matrix in the input *K*. This function implements
% Algorithm 4 from [Górecki et al., 2016a]. Note that the outputted indices
% are ordered according to how much the negative correlation is reduced by
% flipping the corresponding variable, beginning with the one that reduces
% the negative correlation the most.
%
% References:
% [Górecki et al., 2016a] Górecki, J., Hofert, M., and Holeòa, M. (2016). An 
%     approach to structure determination and estimation of hierarchical
%     Archimedean copulas and its application to bayesian classication.
%     Journal of Intelligent Information Systems, pages 21-59.
%
%
% Copyright 2017 Jan Górecki

d = size(K,2);

negatives = K < 0;
if (sum(sum(negatives)) > 0)
    sumCorrOld = -Inf;
    counter = 0;
    
    while (sum(sum(K)) > sumCorrOld)
        sumCorrOld = sum(sum(K));
        newCorr = zeros(1, d);
        
        for i = 1:d
            KNew = K;
            KNew(i, :) = -KNew(i, :);
            KNew(:, i) = -KNew(:, i);
            newCorr(i) = sum(sum(KNew)); %how correlation sum changes after inverting the i-th variable
        end
        invertInd = find(newCorr == max(newCorr), 1, 'first');
        counter = counter + 1;
        invInd(counter) = invertInd;
        
        K(invertInd, :) = -K(invertInd, :);
        K(:, invertInd) = -K(:, invertInd);
    end
else
    invInd = 1;
    counter = 1;
end

out = invInd(1:counter - 1);
