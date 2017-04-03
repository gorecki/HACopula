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
% References:
% [Górecki et al., 2016b] Górecki, J., Hofert, M., and Holeòa, M. (2016). On
%     structure, family and parameter estimation of hierarchical
%     Archimedean copulas. Submitted for publication.
% [Hofert, 2010] Hofert, M. (2010). Sampling Nested Archimedean Copulas
%     with Applications to CDO Pricing. Suedwestdeutscher Verlag fuer
%     Hochschulschriften.
%
%
% Copyright 2017 Jan Górecki

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
        out = @(t, v0x) (1./(theta1 - theta1*log(t + 3060513257434037/1125899906842624).^(1/theta2) + log(t + 3060513257434037/1125899906842624).^(1/theta2)).^v0x);
    case 'C20'
        out = @(t, v0x) (1./exp(v0x.*(1./(1./log(t + 3060513257434037/1125899906842624).^(1/theta2)).^theta1 - 1)));
    case 'FF'
        out = @(t, v0x) ((1-(1-exp(-t).*(1-exp(-theta2))).^(theta1./theta2))/(1-exp(-theta1))).^v0x;
    case 'GG'
        out = @(t, v0x) (exp(-v0x.*t.^(theta1/theta2)));
    case 'C12'
        out = @(t, v0x) (exp(-v0x.* (1./(1./(t.^(1/theta2) + 1)).^theta1 - 1)));
    case 'C14'
        out = @(t, v0x) (exp(-v0x.* (1./(1./(t.^(1./theta2) + 1).^theta2).^theta1 - 1 )));
    case 'JJ'
        out = @(t, v0x) (exp(-v0x.* (-log(1 - ((1 - exp(-t)).^(1./theta2)).^theta1))));
    case '1212'
        out = @(t, v0x) (exp(-v0x.* ((t.^(1./theta2)).^theta1)));
    case '1919'
        out = @(t, v0x) (exp(-v0x.* ((t + exp(theta2)).^(theta1/theta2) - exp(theta1))));   
    otherwise
        error(['get_nonsym_generator01: class combination ' cnames ' not supported']);
end
