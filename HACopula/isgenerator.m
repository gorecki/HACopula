function isGen = isgenerator(family, parameter)
%ISGENERATOR - Check if the combination (family, parameter) corresponds to a
% generator.
%
% If it is a generator, returns 1, otherwise 0. Note that only
% generators with tau >= 0 are concerned, i.e., for example, ('C', -0.5),
% even if it is a generator, it corresponds to tau < 0 and thus the
% function returns 0 - this is connected to the SNC (see [Górecki et
% al., 2016b] for details).
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
        if (parameter >= 0) && (parameter < 1)
            isGen = true;
        else
            isGen = false;
        end
    case {'G', 'J', '12', '14'}
        if (parameter >= 1) && (parameter < Inf)
            isGen = true;
        else
            isGen = false;
        end
    case {'C', 'F', '19', '20'}
        if (parameter > 0) && (parameter < Inf)
            isGen = true;
        else
            isGen = false;
        end
    case '?' % a special arbitrary family that serves only for pre-collapsing
        if (parameter >= -1) && (parameter <= 1)
            isGen = true;
        else
            isGen = false;
        end
    otherwise
        isGen = false;
end
            
    