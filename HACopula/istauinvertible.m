function invertible = istauinvertible(family, tauMat)
%ISTAUINVERTIBLE - returns 1, if all taus from the matrix of taus *tauMat*
%could be inverted through the theta-tau relationship assuming the family
%*family*. Otherwise, returns 0.
%
% References:
% [Górecki et al., 2016b] Górecki, J., Hofert, M., and Holeòa, M. (2016). On
%     structure, family and parameter estimation of hierarchical
%     Archimedean copulas. arXiv preprint arXiv:1611.09225.
%
%
% Copyright 2017 Jan Górecki

invertible = zeros(size(tauMat));
for i = 1:size(tauMat,1)
    for j = 1:size(tauMat,2)
        invertible(i, j) = istauinvertiblescalar(family, tauMat(i,j));
    end
end
invertible = prod(reshape(invertible,[1, numel(tauMat)]));


function invertible = istauinvertiblescalar(family, tau)
% returns true iif *tau* lies in the interval \tau_{(a)}(\Theta_a), see [Górecki et al., 2016b]

% \tau_{(a)}(\Theta_a)
tauRange = getfamilytaurange(family);


% true iif theta lies in the interval tauTheta
invertible = min(max(tauRange(1),tau),tauRange(2)) == tau;
        
        
        
        