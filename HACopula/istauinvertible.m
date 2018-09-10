function invertible = istauinvertible(family, tauMat)
%ISTAUINVERTIBLE - returns 1, if all taus from the matrix of taus *tauMat*
%could be inverted through the theta-tau relationship assuming the family
%*family*. Otherwise, returns 0.
%
% References:
% [Gorecki et al., 2017] On Structure, Family and Parameter Estimation
%     of Hierarchical Archimedean copulas. Journal of Statistical Computation 
%     and Simulation, 87(17), 3261-3324
%
%
% Copyright 2018 Jan Gorecki

invertible = zeros(size(tauMat));
for i = 1:size(tauMat,1)
    for j = 1:size(tauMat,2)
        invertible(i, j) = istauinvertiblescalar(family, tauMat(i,j));
    end
end
invertible = prod(reshape(invertible,[1, numel(tauMat)]));

end


function invertible = istauinvertiblescalar(family, tau)
% returns true iif *tau* lies in the interval \tau_{(a)}(\Theta_a), see [Gorecki et al., 2017]

% \tau_{(a)}(\Theta_a)
tauRange = getfamilytaurange(family);


% true iif theta lies in the interval tauTheta
invertible = min(max(tauRange(1),tau),tauRange(2)) == tau;
        
end        
        
        
