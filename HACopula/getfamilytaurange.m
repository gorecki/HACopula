function familyRange = getfamilytaurange(family)
%GETFAMILYTAURANGE - returns the parameter range of a given family
%
% References:
% [Górecki et al., 2016b] Górecki, J., Hofert, M., and Holeòa, M. (2016). On
%     structure, family and parameter estimation of hierarchical
%     Archimedean copulas. Submitted for publication.
%
%
% Copyright 2017 Jan Górecki

switch family
    case 'A'
        familyRange = [0 1/3-eps(1/3)]; % [0, 1/3)
    case {'C', 'F', '20'}
        familyRange = [eps(0) 1-eps(1)]; % (0, 1)
    case {'G', 'J'}
        familyRange = [0 1-eps(1)]; % [0, 1)
    case {'12', '14'}
        % here, it should be [1/3, 1), however, due to limited precision of
        % double, the lower bound is shifted to 1/3+eps(1/3) because 1/3 in
        % Matlab is lower than real (infinitely precise) 1/3
        familyRange = [1/3+eps(1/3) 1-eps(1)];  
    case '19'
        familyRange = [1/3+eps(1/3) 1-eps(1)]; % (1/3, 1)
    case '?'
        familyRange = [-1 1]; % an arbitrary interval that corresponds to the maximal range of Kendall's tau          
    otherwise
        error('istauinvertible: Unsupported family.')
end