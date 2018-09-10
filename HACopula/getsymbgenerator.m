function out = getsymbgenerator(family, inversion)
%GETSYMBGENERATOR - returns a symbolic version of an Archimedean generator
%function.
%
% out = getsymbgenerator(family, inverted) returns the function of the
% Archimedean generator from the input family. If
% inverted == 1, then the inverse (\psi^{-1}) of the function is returned,
% otherwise, the original function \psi is returned.
% 
% NOTE:
% In Octave, the OctSymPy package is needed (https://github.com/cbm755/octsympy).  
%
% References:
% [Nelsen, 2006] Nelsen, R. (2006). An Introduction to Copulas. Springer,
% 2nd edition.
%
%
% Copyright 2018 Jan Gorecki

syms t;
syms theta;

if (inversion == 0)
    %psi
    switch family
        case 'A'
            out = (1 - theta)/(exp(t) - theta);
        case 'C'
            out = (1 + t) ^ (-1/theta);
        case 'F'
            out = -log(1 - (1 - exp(-theta)) * exp(-t))/theta;
        case 'G'
            out = exp(-t ^ (1/theta));
        case 'J'
            out = 1 - (1 - exp(-t)) ^ (1/theta);
        case '12'
            out = 1/(1 + t^(1/theta));
        case '14'
            out = (1 + t^(1/theta))^(-theta);
        case '19'
            out = theta/log(t + exp(theta));
        case '20'
            out = (log(t + exp(1)))^(-1/theta);
        otherwise
            error('HACopula:BadInputs', 'getsymbgenerator: Unknown Archimedean family.');
    end
else
    %the inversion of the generator
    switch family
        case 'A'
            out = log((1 - theta  * (1 - t)) / t);
        case 'C'
            out = t ^ (-theta) - 1;
        case 'F'
            out = -log ((1 - exp(-theta * t))/(1 - exp(-theta)));
        case 'G'
            out = (-log(t))^theta;
        case 'J'
            out = -log(1 - (1 - t)^theta);
        case '12'
            out = (1/t - 1)^theta;
        case '14'
            out = (t^(-1/theta) - 1)^theta;
        case '19'
            out = exp(theta/t) - exp(theta);
        case '20'
            out = exp(t^(-theta)) - exp(1);
        otherwise
            error('HACopula:BadInputs', 'getsymbgenerator: Unknown Archimedean family');
    end
end
    
end
