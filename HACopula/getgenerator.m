function out = getgenerator(family, theta, inverted)
%GETGENERATOR - returns an Archimedean generator function
%
% out = getgenerator(family, theta, inverted) returns the function of the
% Archimedean generator from the input family with the parameter theta. If
% inverted == 1, then the inverse (\psi^{-1}) of the function is returned,
% otherwise, the original function \psi is returned.
%
% References:
% [Nelsen, 2006] Nelsen, R. (2006). An Introduction to Copulas. Springer,
% 2nd edition.
%
%
% Copyright 2017 Jan Górecki



if (inverted == 0)
    % psi
    switch family
        case 'A'
            out = @(t)((1 - theta)./(exp(t) - theta));
        case 'C'
            out = @(t)((1 + t) .^ (-1/theta));
        case 'F'
            out = @(t)(-log(1 - (1 - exp(-theta)).*exp(-t))./theta);
        case 'G'
            out = @(t)(exp(-t.^(1/theta)));
        case 'J'
            out = @(t)(1 - (1 - exp(-t)).^(1/theta));
        case '12'
            out = @(t)(1./(1 + t.^(1/theta)));
        case '14'
            out = @(t)((1 + t.^(1/theta)).^(-theta));
        case '19'
            out = @(t)(theta./log(t + exp(theta)));
        case '20'
            out = @(t)(log(t + exp(1)).^(-1/theta));
        otherwise
            error(['getgenerator: The family ''' family ''' is not supported.']);
    end
else
    % psi^{-1}
    switch family
        case 'A'
            out = @(t)(log((1 - theta  * (1 - t)) ./ t));
        case 'C'
            out = @(t)(t .^ (-theta) - 1);
        case 'F'
            out = @(t)(-log ((1 - exp(-theta .* t))./(1 - exp(-theta))));
        case 'G'
            out = @(t)((-log(t)).^theta);
        case 'J'
            out = @(t)(log(1./ ( 1 - (1 - t).^theta)));
        case '12'
            out = @(t)((1./t - 1).^theta);
        case '14'
            out = @(t)((t.^(-1/theta) - 1).^theta);
        case '19'
            out = @(t)(exp(theta./t) - exp(theta));
        case '20'
            out = @(t)(exp(t.^(-theta)) - exp(1));
        otherwise
            error(['getgenerator: The family ''' family ''' is not supported.']);
    end
end