function coefficient = gettaildependence(family, theta, type)
% GETTAILDEPENDENCE - upper- and lower-tail dependence coefficitents.
% 
% Purpose:
% Computes the upper- and lower-tail dependence coefficitents according to [Nelsen, 2006]. 
%
% Inputs:
% family    - a family of Archimedean generator
% theta     - the paramter of the generator
% type      - 'lower' or 'upper' for the correspoding type of tail
%             coefficient
%
% Output:
% coefficient - the value of the corresponding coefficient
%
% References:
% [Nelsen, 2006] Nelsen, R. (2006). An Introduction to Copulas. Springer,
%     2nd edition.
%
% Copyright 2018 Jan Gorecki

if strcmp(type,'lower')
    switch family
        case {'A', 'F', 'G', 'J'}
            coefficient = 0;
        case {'C', '12'}
            coefficient = 2^(-1/theta);
        case '14'
            coefficient = 0.5;
        case {'19', '20'}
            coefficient = 1;
        otherwise
            error('HACopula:BadInputs', 'gettaildependence: Unknown family.');
    end
elseif strcmp(type,'upper')
    switch family
        case {'A', 'C', 'F', '19', '20'}
            coefficient = 0;
        case {'G', 'J', '12', '14'}
            coefficient = 2 - 2^(1/theta);
        otherwise
            error('HACopula:BadInputs', 'gettaildependence: Unknown family.');
    end
    
else
    error('HACopula:BadInputs', 'gettaildependence: Unsupported type of a tail. Choose from {''lower'', ''upper''}');
end

end