function isGen = isgenerator(family, parameter)
%ISGENERATOR - Check if the combination (family, parameter) corresponds to a
% generator.
%
% If it is a generator, returns 1, otherwise 0. Note that only
% generators with tau >= 0 are concerned, i.e., for example, ('C', -0.5),
% even if it is a generator, it corresponds to tau < 0 and thus the
% function returns 0 - this is connected to the SNC (see [Gorecki et al., 2017] 
% for details).
%
% References:
% [Gorecki et al., 2017] On Structure, Family and Parameter Estimation
%     of Hierarchical Archimedean copulas. Journal of Statistical Computation 
%     and Simulation, 87(17), 3261ÿ3324
%
%
% Copyright 2018 Jan Gorecki

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
   
end
    