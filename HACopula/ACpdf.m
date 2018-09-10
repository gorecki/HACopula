function pdf = ACpdf(family, theta, u1, u2)
% ACPDF - evaluates the density function of a bivariate Archimedean copula
%
% pdf = ACpdf(family, theta, u1, u2) evaluates the pdf of the 2-AC from FAMILY
% with parameter THETA at points (u1, u2)
%
% Inputs:
% family    - 'A', 'C', 'F', 'G', 'J', '12', '14', '19' or '20'
% theta     - the parameter of the copula family
% u1, u2    - two [0, 1]-valued vectors of the same length 
%
% Output:
% pdf       - the density evaluated at (u1, u2)
%
% References:
% [Joe, 2014] Joe, H. (2014) Dependence modeling with copulas. CRC Press.
%
% Copyright 2018 Jan Gorecki

% NOTE:
% When introducing a new family (new density function), e.g., from
% literature, it is very helpful to first generate the density formula
% using, e.g.:
% 0) syms u1; syms u2; syms theta3;
% 1) AC = HACopula({{'20', 1.5}, 1, 2});
% 2) simplify(diff(diff(getcdf(AC), u1), u2))
%
% and then check if this function gives the same values as the one taken
% from the literature.
%
% All formulas below are taken from [Joe, 2014].


switch family
    case 'A'
        u12 = theta .* (1 - u1) .* (1 - u2);
        pdf = (1 - u12).^(-3) .* (1 - theta + 2 .* theta .* u1 .* u2 - (1 - theta) .* u12);
    case 'C'
        pdf = (1 + theta) .* (u1 .* u2).^(-theta-1) .* (u1.^(-theta) + u2.^(-theta) - 1).^(-2-1/theta);
    case 'F'
        pdf = (-theta .* expm1(-theta) .* exp(-theta .* (u1 + u2))) ./ (-expm1(-theta) - expm1(-theta .* u1) .* expm1(-theta .* u2) ).^2;
    case 'G'
        x = -log(u1);
        y = -log(u2);
        xy = x.^theta + y.^theta;
        pdf = exp(-xy.^(1/theta)) .* (xy.^(1/theta) + theta - 1) .* xy.^(1/theta-2) .* (x .* y).^(theta-1) ./ (u1 .* u2);
    case 'J'
        u1t = (1 - u1).^theta;
        u2t = (1 - u2).^theta;
        u12t = u1t + u2t - u1t .* u2t;
        pdf = u12t.^(1/theta-2) .* (1 - u1).^(theta-1) .* (1 - u2).^(theta-1) .* (theta - 1 + u12t);
    case '12'
        pdf = BB1(u1, u2, 1, theta);
    case '14'
        pdf = BB1(u1, u2, 1/theta, theta);
    case '19'
        pdf = BB2(u1, u2, 1, theta);
    case '20'
        pdf = BB2(u1, u2, theta, 1);
    otherwise
        error('HACopula:ACpdf', ['ACpdf: Family ''' family ''' is not supported.']);
end

end

function pdf = BB1(u1, u2, theta, delta)
x = (u1.^(-theta) - 1).^delta;
y = (u2.^(-theta) - 1).^delta;

pdf = (1+(x+y).^(1/delta)).^(-1/theta-2) .* (x+y).^(1/delta-2) .* ...
      (theta * (delta - 1) + (theta * delta + 1) .* (x + y).^(1/delta)) .* ...
      (x .* y).^(1-1/delta) .* (u1 .* u2).^(-theta-1);

end

function pdf = BB2(u1, u2, theta, delta)
x = expm1(delta .* (u1.^(-theta) - 1));
y = expm1(delta .* (u2.^(-theta) - 1));

pdf = (1 + 1/delta .* log(x + y + 1)).^(-1/theta-2) .* (x + y + 1).^(-2) .* ...
      (1 + theta + theta * delta + theta .* log(x + y + 1)) .* ...
      (x + 1) .* (y + 1) .* (u1 .* u2).^(-theta-1);
end
