function thetaEst = mlefor2AC(U, family, tau, lBound, uBound)
%MLEFOR2AC - Maximum likelihood estimation for a 2-AC
%
% thetaEst = mlefor2AC(U, family, tau, lBound, uBound) returns the maximum
% likelihood estimate for the 2-AC from the input family based on the
% observations U.
%
% Inputs:
% U         - An (n * 2)-matrix with observations from [0,1]^2.
% family    - 'A', 'C', 'F', 'G', 'J', '12', '14', '19' or '20'.
% tau       - The estimate of Kendall's based on U, which is used as an
%             initial value (transformed to the corresponding parameter
%             value) for MLE. If tau does not relate to any parameter from
%             the given family, (lBound + uBound)/2 is used as the initial
%             value provided both lBound and uBound are not Inf. If lBound
%             or uBound is Inf, then the initial value is set to 1 - SHIFT,
%             where SHIFT is defined below.
% lBound    - Lower bound for the output parameter.
% uBound    - Upper bound for the output parameter.
%
% NOTE: If a NaN or an Inf is generated in the estimation process, it is
% approximated using the nanapprox function.
%
% References:
% [Genest and Favre, 2007] Genest, C. and Favre, A. (2007). Everything you always
% wanted to know about copula modeling but were afraid to ask. Hydrol. Eng.,
% 12:347-368.
%
%
% Copyright 2017 Jan Górecki

% The densities below for the given family are obtained through:
% 0) syms u1; syms u2; syms theta3;
% 1) AC = HACopula({{'20', 1.5}, 1, 2});
% 2) simple(diff(diff(getcdf(AC), u1), u2))

switch family
    case 'A'
        pdf = @(u1, u2, theta) ((u2.*(u1 - 1) - u1 + 1).*theta.^2 + (u1 + u2.*(u1 + 1) - 2).*theta + 1)./(theta.*u1 - theta + theta.*u2 - theta.*u1.*u2 + 1).^3;
    case 'C'
        pdf = @(u1, u2, theta) (theta + 1)./(u1.^(theta + 1).*u2.^(theta + 1).*(1./u1.^theta + 1./u2.^theta - 1).^(1./theta + 2));
    case 'F'
        pdf = @(u1, u2, theta) -(theta.*(exp(theta.*(u1 + u2 + 1)) - exp(theta.*(u1 + u2 + 2))))./(exp(2.*theta.*(u1 + 1)) - ...
            2.*exp(theta.*(u2 + 2)) - 2.*exp(theta.*(u1 + 2)) + exp(2.*theta.*(u2 + 1)) + exp(2.*theta) + 2.*exp(theta.*(u1 + u2 + 1))...
            + 2.*exp(theta.*(u1 + u2 + 2)) + exp(2.*theta.*(u1 + u2)) - 2.*exp(theta.*(u1 + 2.*u2 + 1)) - 2.*exp(theta.*(2.*u1 + u2 + 1)));
    case 'G'
        pdf = @(u1, u2, theta) (exp(-((-log(u1)).^theta + (-log(u2)).^theta).^(1./theta)).*(-log(u1)).^theta.*(-log(u2)).^theta.*((-log(u1)).^theta +...
            (-log(u2)).^theta).^(1./theta).*(theta + ((-log(u1)).^theta + (-log(u2)).^theta).^(1./theta) - 1))./...
            (u1.*u2.*log(u1).*log(u2).*((-log(u1)).^theta + (-log(u2)).^theta).^2);
    case 'J'
        pdf = @(u1, u2, theta) ((1 - u1).^theta.*(1 - u2).^theta.*((1 - u1).^theta + (1 - u2).^theta - (1 - u1).^theta.*(1 - u2).^theta).^(1./theta - 2).*...
            (theta + (1 - u1).^theta + (1 - u2).^theta - (1 - u1).^theta.*(1 - u2).^theta - 1))./((u1 - 1).*(u2 - 1));
    case '12'
        pdf = @(u1, u2, theta) ((1./u1 - 1).^theta.*(1./u2 - 1).^theta.*(((1./u1 - 1).^theta + (1./u2 - 1).^theta).^(2./theta).*(theta + 1) +...
            ((1./u1 - 1).^theta + (1./u2 - 1).^theta).^(1./theta).*(theta - 1)))./(u1.*u2.*(((1./u1 - 1).^theta + (1./u2 - 1).^theta).^...
            (1./theta) + 1).^3.*((1./u1 - 1).^theta + (1./u2 - 1).^theta).^2.*(u1 - 1).*(u2 - 1));
    case '14'
        pdf = @(u1, u2, theta) ((1./u1.^(1./theta) - 1).^theta.*(1./u2.^(1./theta) - 1).^theta.*((1./u1.^(1./theta) - 1).^theta + (1./u2.^(1./theta) - 1).^theta).^...
            (1./theta - 2).*(theta + 2.*theta.*((1./u1.^(1./theta) - 1).^theta + (1./u2.^(1./theta) - 1).^theta).^(1./theta) - 1))./...
            (theta.*u1.*u2.*(u1.^(1./theta) - 1).*(u2.^(1./theta) - 1).*(((1./u1.^(1./theta) - 1).^theta + (1./u2.^(1./theta) - 1).^theta).^(1./theta) + 1).^(theta + 2));
    case '19'
        pdf = @(u1, u2, theta) (theta.^3.*exp(theta./u1 + theta./u2).*(log(exp(theta./u1) + exp(theta./u2) - exp(theta)) + 2))./...
            (u1.^2.*u2.^2.*log(exp(theta./u1) + exp(theta./u2) - exp(theta)).^3.*(exp(theta./u1) + exp(theta./u2) - exp(theta)).^2);
    case '20'
        pdf = @(u1, u2, theta) (1267650600228229401496703205376.*exp(1./u1.^theta + 1./u2.^theta).*(theta + theta.*log(exp(1./u1.^theta) + exp(1./u2.^theta)...
            - 3060513257434037./1125899906842624) + 1))./(u1.^(theta + 1).*u2.^(theta + 1).*log(exp(1./u1.^theta) + exp(1./u2.^theta)...
            - 3060513257434037./1125899906842624).^(1./theta + 2).*(1125899906842624.*exp(1./u1.^theta) + 1125899906842624.*exp(1./u2.^theta) - 3060513257434037).^2);
    otherwise
        error(['mlefor2AC: Family ''' family ''' is not supported.']);
end

if istauinvertible(family, tau) 
    % is tau inversion estimate if possible
    start = tau2theta(family, tau);
elseif ~isinf(lBound) && ~isinf(uBound)
    % use the mean of the parameter range
    start = mean([lBound uBound]);  % start in the midle of [lBound, uBound]
else
    % use predefined starting values 
    SHIFT = 0.1;
    switch family
        case {'G', 'J', '12', '14'}
            start = 1 + SHIFT;
        case {'C', 'F', '19', '20'}
            start = SHIFT;
        case 'A'
            start = 1 - SHIFT;
    end
end


FMLE = @(theta)(-sum(log(pdfNaNCheck(pdf(U(:,1), U(:,2), max(lBound,min(uBound,theta))), U)))); % do NaN check and also bound (trim) theta to [lBound, uBound]
opts = optimset('MaxFunEvals', 20, 'Display', 'off');
%opts = optimset('Display', 'off');
%thetaEst = fminbnd(FMLE, lBound, min(uBound, 1000), opts);
thetaEst = fminsearch(FMLE, start, opts);
end

function likelihoods = pdfNaNCheck(likelihoods, data)
% NaN check
[likelihoods, nNaNs] = nanapprox(likelihoods, data);
if nNaNs > 0
    if nNaNs < length(likelihoods)
        warning(['mlefor2AC::fminsearch::nanapprox: ' num2str(nNaNs) ' NaNs detected. Replacing by the likelihoods of the closest (in the Euclidian distance) data with non-NaN likelihoods.']);
    else
        warning(['mlefor2AC::fminsearch::nanapprox: ' num2str(nNaNs) ' NaNs detected. Trying another theta value.']);
    end
end
end

