function out = getnestedgenerators(family1, family2, theta1, theta2)
%GETNESTEDGENERATORS - Gets an inner generator.
%
% out = getnestedgenerators(family1, family2, theta1, theta2) returns the 
% inner generators, i.e., exp(-V_0 * \psi_0^{-1}(\psi_1(t))), where \psi_0
% is from family1 with the parameter theta 1 and where \psi_1
% is from family2 with the parameter theta 2.
%
% NOTE: The functions in the code below are generated using:
% syms theta1;
% syms theta2;
% syms t;
% syms theta;
% %\psi_0^{-1}(\psi_1(t)) = 
% simplify(subs(subs(getsymbgenerator('19',1), theta, theta1), t, subs(getsymbgenerator('19',0), theta, theta2)))
%
% NOTE2:
% In some of the formulas, exp(x)-1 is substituted by more accurate expm1(x)
%
%
% References:
% [Gorecki et al., 2017] On Structure, Family and Parameter Estimation
%     of Hierarchical Archimedean copulas. Journal of Statistical Computation 
%     and Simulation, 87(17), 3261ÿ3324
% [Hofert, 2010] Hofert, M. (2010). Sampling Nested Archimedean Copulas
%     with Applications to CDO Pricing. Suedwestdeutscher Verlag fuer
%     Hochschulschriften.
%
%
% Copyright 2018 Jan Gorecki

cnames = [family1 family2];

switch cnames
    case 'AA'
        out = @(t, v0x) exp(-v0x .*(log(-(theta1 - theta2 + exp(t) - theta1.*exp(t))./(theta2 - 1))));
    case 'AC'
        out = @(t, v0x)(  (1 + (1 - theta1) * ( (1 + t) .^ (1/theta2) - 1 )).^(-v0x)  );
    case 'A19'
        out = @(t, v0x) ( (theta1 - (log(t + exp(theta2))*(theta1 - 1))/theta2) .^(-v0x) );
    case 'C19'
        out = @(t, v0x) (  exp(-v0x .* (1./(theta2./log(t + exp(theta2))).^theta1 - 1) ) );
    case 'CC'
        out = @(t, v0x) ( exp( -v0x .* ((1+t).^(theta1/theta2) - 1) ) );
    case 'A20'
        out = @(t, v0x) (1./(theta1 - theta1*log(t + exp(1)).^(1/theta2) + log(t + exp(1)).^(1/theta2)).^v0x);
    case 'C20'
        out = @(t, v0x) (1./exp(v0x.*(1./(1./log(t + exp(1)).^(1/theta2)).^theta1 - 1)));
    case 'FF'
        out = @(t, v0x) ((1-(1-exp(-t).*(-expm1(-theta2))).^(theta1./theta2))/(-expm1(-theta1))).^v0x;
    case 'GG'
        out = @(t, v0x) (exp(-v0x.*t.^(theta1/theta2)));
    case 'C12'
        out = @(t, v0x) (exp(-v0x.* (1./(1./(t.^(1/theta2) + 1)).^theta1 - 1)));
    case 'C14'
        out = @(t, v0x) (exp(-v0x.* (1./(1./(t.^(1./theta2) + 1).^theta2).^theta1 - 1 )));
    case 'JJ'
        out = @(t, v0x) (exp(-v0x.* (-log(1 - ((-expm1(-t)).^(1./theta2)).^theta1))));
    case '1212'
        out = @(t, v0x) (exp(-v0x.* ((t.^(1./theta2)).^theta1)));
    case '1919'
        out = @(t, v0x) (exp(-v0x.* ((t + exp(theta2)).^(theta1/theta2) - exp(theta1))));   
    otherwise
        error('HACopula:getnestedgenerators', ['getnestedgenerators: class combination ' cnames ' not supported']);
end

end