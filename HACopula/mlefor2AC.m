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
% Copyright 2018 Jan Gorecki

% The densities below for the given family are obtained through:
% 0) syms u1; syms u2; syms theta3;
% 1) AC = HACopula({{'20', 1.5}, 1, 2});
% 2) simplify(diff(diff(getcdf(AC), u1), u2))


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

pdf = @(u1, u2, theta) ACpdf(family, theta, u1, u2);
FMLE = @(theta)(-sum(log(pdfNaNCheck(pdf(U(:,1), U(:,2), max(lBound,min(uBound,theta))), U)))); % do NaN check and also bound (trim) theta to [lBound, uBound]
%opts = optimset('MaxFunEvals', 20, 'Display', 'off');
opts = optimset('Display', 'off');
%thetaEst = fminbnd(FMLE, lBound, min(uBound, 1000), opts);
thetaEst = fminsearch(FMLE, start, opts);
end

function likelihoods = pdfNaNCheck(likelihoods, data)
% NaN check
[likelihoods, nNaNs] = nanapprox(likelihoods, data);
if nNaNs > 0
    if nNaNs < length(likelihoods)
        warning('HACopula:NaN_detected', ['mlefor2AC::fminsearch::nanapprox: ' num2str(nNaNs) ' NaNs detected. Replacing by the likelihoods of the closest (in the Euclidian distance) data with non-NaN likelihoods.']);
    else
        warning('HACopula:mlefor2AC', ['mlefor2AC::fminsearch::nanapprox: ' num2str(nNaNs) ' NaNs detected. Trying another theta value.']);
    end
end
end

