function cCdf = ACcondcdf(family, theta, u1, u2)
%CONDCDF - evalueates the conditional CDF for a bivariate Archimedean copula
%
% cCdf = ACcondcdf(family, theta, u1, u2) evaluates the conditional CDF 
% (\delta C / \delta u1, commonly denoted by C(u2|u1)) of the 2-AC from FAMILY
% with parameter THETA at points (u1, u2)
%
% Inputs:
% family    - 'A', 'C', 'F', 'G', 'J', '12', '14', '19' or '20'
% theta     - the parameter of the copula family
% u1, u2    - two [0, 1]-valued vectors of the same length 
%
% Output:
% cCdf       - the conditional density evaluated at (u1, u2)
%
% References:
% [Joe, 2014] Joe, H. (2014) Dependence modeling with copulas. CRC Press.
%
% Copyright 2018 Jan Gorecki

% NOTE:
% When introducing a new family (new density function), e.g., from
% literature, it is very helpful to first generate the conditional CDF formula
% using, e.g.:
% 0) syms u1; syms u2; syms theta3;
% 1) AC = HACopula({{'20', 1.5}, 1, 2});
% 2) simplify(diff(getcdf(AC), u1))
%
% and then check if this function gives the same values (outputs) as the one taken
% from the literature.
%
% All formulas below are taken from [Joe, 2014].


 switch family
        case 'A'
            cCdf = BB10(u1, u2, 1, theta);
        case 'C'
            cCdf = (1 + u1.^(theta) .* (1./u2.^theta - 1)).^(-1-1/theta);
        case 'F'
            cCdf = exp(-theta .* u1) ./ (expm1(-theta) ./ expm1(-theta .* u2) + expm1(-theta .* u1));
        case 'G'
            x = -log(u1);
            y = -log(u2);
            cCdf = exp(-(x.^theta + y.^theta).^(1/theta)) .* (1+(y./x).^theta).^(1/theta-1) ./ u1;
        case 'J'
            cCdf = (1 + ((1 - u2) ./ (1 - u1)).^theta - (1 - u2).^theta).^(-1+1/theta) .* (1 - (1 - u2).^theta);
        case '12'
            cCdf = BB1(u1, u2, 1, theta);
        case '14'
            cCdf = BB1(u1, u2, 1/theta, theta);
        case '19'
            cCdf = BB2(u1, u2, 1, theta);
        case '20'
            cCdf = BB2(u1, u2, theta, 1);
        otherwise
            error('HACopula:ACcondcdf', 'ACcondcdf: Unsupported family.');
 end
    
end

function cCdf = BB1(u1, u2, theta, delta)
x = (u1.^(-theta) - 1).^delta;
y = (u2.^(-theta) - 1).^delta;

cCdf = (1+(x+y).^(1/delta)).^(-1/theta-1) .* (x+y).^(1/delta-1) .* x.^(1-1/delta) .* u1.^(-theta-1);
end


function cCdf = BB2(u1, u2, theta, delta)
x = expm1(delta .* (u1.^(-theta) - 1));
y = expm1(delta .* (u2.^(-theta) - 1));

cCdf = (1 + 1/delta .* log(x+y+1)).^(-1/theta-1) ./ (x+y+1) .* (x+1) .* u1.^(-theta-1);
end

function cCdf = BB10(u1, u2, theta, delta)
cCdf = (1 - delta .* (1 - u1.^theta) .* (1 - u2.^theta)).^(-1/theta-1) .* u2 .* (1 - delta .* (1 - u2.^theta));
end

