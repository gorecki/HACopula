function kendTrans = ACkendtrans(family, theta, t)
%ACKENDTRANS - evaluates the Kendall transformation of a bivariate Archimedean copula
%
% kendTrans = ACkendtrans(family, theta, t) evaluates the Bivariate
% probability integral transform (known as the Kendall transformation) of
% the 2-AC C from FAMILY with parameter THETA at input T, i.e., given that
% W = C(U1,U2) and assuming that K is the distribution function of W, then
% kendTrans = K(t).
%
% Inputs:
% family    - 'A', 'C', 'F', 'G', 'J', '12', '14', '19' or '20'
% theta     - the parameter of the copula family
% t         - a [0, 1]-valued vector of inputs for the evaluation
%
% Output:
% kendTrans - the the Bivariate probability integral transform evaluated at t
%
% References:
% [Genest and Favre, 2007] Genest, C. and Favre, A. (2007). Everything you always
% wanted to know about copula modeling but were afraid to ask. Hydrol. Eng.,
% 12:347-368.
%
%
% Copyright 2018 Jan Gorecki

% NOTE:
% the formulas below are basically computed using:
% 0) syms t; syms theta; family = 'A';
% 1) simplify(t - getsymbgenerator(family, 1) / diff(getsymbgenerator(family, 1), t))
%
% In some of the formulas, exp(x)-1 is substituted by more accurate
% expm1(x) and some repeating terms are extracted and computed beforehand


switch family
    case 'A'
        tTh = t.*theta - theta + 1;
        kendTrans = t - (t .* log(tTh ./ t) .* tTh)./(theta - 1);
    case 'C'
        kendTrans = t - (t .* (t.^theta - 1)) ./ theta;
    case 'F'
        em1T = expm1(-t .* theta);
        kendTrans = t + (exp(t .* theta) .* log(em1T ./ expm1(-theta)) .* em1T) ./ theta;
    case 'G'
        kendTrans = t - (t.*log(t))./theta;
    case 'J'
        tT = (1 - t).^theta;
        kendTrans = t + (log(1 - tT) .* (tT - 1) .* (1 - t).^(1 - theta)) ./ theta;
    case '12'
        kendTrans = t - (t.*(t - 1))./theta;
    case '14'
        kendTrans = 2 .* t - t.^(1/theta + 1);
    case '19'
        kendTrans = t - (t.^2 .* (expm1((theta .* (t - 1)) ./t ))) ./ theta;
    case '20'
        tT = t.^(theta + 1);
        kendTrans = t - ((exp(1) .* tT .* exp(-1./t.^theta)) - tT)./theta;
    otherwise
        error('HACopula:ACkendtrans', 'ACkendtrans: Unsupported family.');
end

end