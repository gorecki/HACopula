function [areEqual, ratio] = comparestructures(HAC1, HAC2)
%COMPARESTRUCTURES - compare two HAC structures 
%
% areEqual = comparestructures(HAC1, HAC2) returns 1, if the structures of
% HAC1 and HAC2 are equal, otherwise, returns 0.
%
% [areEqual, ratio] = comparestructures(HAC1, HAC2) additionaly returns a
% ratio expressing a certain type of similarity between the two structures.
% This ratio is based on the trivariate decompositon of a HAC tree
% structure introduced in [Segers and Uyttendaele (2014)]
%
% NOTE:
% If one is just interested if the two structures are equal (and not in the
% ratio), one should require only the first output (areEqual), for which
% the computation method is quite fast. Otherwise, i.e., if the ratio is
% also required, a more computionally demanding method (based on the
% trivariate decompostion) is used.
%
% References:
% [Segers and Uyttendaele (2014)] Segers J, Uyttendaele N (2014).
% "Nonparametric Estimation of the Tree Structure of a Nested Archimedean
% Copula." Computational Statistics & Data Analysis, 72, 190–204. doi:
% 10.1016/j.csda.2013.10.028.
%
%
% Copyright 2018 Jan Gorecki

% check for the same dimensions
if (HAC1.Dim ~= HAC2.Dim)
    areEqual = false;
    ratio = NaN;
    return;
end

if nargout ~= 2 
    %% fast method - just checks the equality, and does not compute the ratio
    
    % here, the set of leaves is taken for each fork of the first HAC, and
    % the resulting set of those sets is checked to match with the
    % corresponding set built for the second HAC.
        
    d = HAC1.Dim;
    
    forks1 = HAC1.Forks;
    forks2 = HAC2.Forks;
    
    if (size(forks1,2) ~= size(forks2,2))
        areEqual = false;
        return;
    end
    
    k = size(forks1,2);
    
    % get leaves of each of the forks1
    leaves1 = cell(1, k);
    leaves1Extend = zeros(k, d);
    leaves2 = cell(1, k);
    leaves2Extend = zeros(k, d);
    for i = 1:size(forks1,2)
        leaves1{i} = forks1{i}.Leaves;
        leaves2{i} = forks2{i}.Leaves;
        
        % extend leaves for sortrow below
        leaves1Extend(i, 1:d) =  [zeros(1, d - size(leaves1{i},2)) sort(leaves1{i})];
        leaves2Extend(i, 1:d) =  [zeros(1, d - size(leaves2{i},2)) sort(leaves2{i})];
    end
    
    leaves1Sorted = sortrows(leaves1Extend);
    leaves2Sorted = sortrows(leaves2Extend);
    areEqual = (sum(sum(leaves1Sorted == leaves2Sorted)) == k*d);

else % nargout == 2 
    %% trivariate approach - slower but gives the ratio

    % here, for each possible triplet (nchoosek(d, 3)) of the leaves, the
    % trivariate HAC in the sense of [Segers and Uyttendaele (2014)] is
    % built/analyzed for each of the two input HACs and these are then checked to
    % match.
    %
    % Given a triplet (i, j, k), according to [Segers and Uyttendaele
    % (2014)], there can exist only four differently structured trivariate
    % HACs. The particular one can be detected using the "youngest common
    % antecendent" (yca (lca in [Segers and Uyttendaele (2014)])).
    % For the pairs (i, j), (i, k) and (j, k), the corresponding ycas can be:
    % 1) (AAB)  2) (ABA)  3) (BAA)  4) (AAA),
    % where A and B belong to the set of the TauOrdering of all forks.
    % To distinguish among these four cases, just two bits are needed:
    % I) yca(i, j) == yca(i, k) and II) yca(i, k) == yca(j, k)
    % The four cases are then identified:
    % 1) = 10  (first bit is from I) and the second from II))
    % 2) = 00
    % 3) = 01
    % 4) = 11
    %
    % Hence, to check the match of two trivariate structures, it is enough
    % to compare these two bits obtained for the first HAC to the
    % corresponding bits obtained for the second HAC. 
    
    nEqual = 0;
    d = HAC1.Dim;
    ycaMat1 = ycatauorderingmatrix(HAC1);
    ycaMat2 = ycatauorderingmatrix(HAC2);
    for i = 1:d
        for j = i+1:d
            for k = j+1:d
                ij_ik1 = ycaMat1(i, j) == ycaMat1(i, k);
                ik_jk1 = ycaMat1(i, k) == ycaMat1(j, k);
                ij_ik2 = ycaMat2(i, j) == ycaMat2(i, k);
                ik_jk2 = ycaMat2(i, k) == ycaMat2(j, k);
                areEqual = all([ij_ik1 ik_jk1] == [ij_ik2 ik_jk2]);
                nEqual = nEqual + areEqual;
            end
        end
    end
    
    ratio = nEqual/nchoosek(d, 3);
    areEqual = ratio == 1;
end

end
