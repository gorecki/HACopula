function Sn = gof2Sn(U, family, theta, gofTestName, emp2copula)
%GOF2SN - Computing a goodness-of-fit statistic for a 2-AC.
%
% Sn = gof2Sn(data, family, parameter, gofTestName, emp2copula) returns the
% goodness-of-fit statistic corresponding to gofTestName ('E', 'K' or 'R') computed for the
% 2-AC from the input family and with the parameter theta and the bivariate
% observations U. To improve efficiency of computation for gofTestName from
% {'E', 'K'}, provide emp2copula computed by computeallemp2copulas(U).
%
% NOTE: if a NaN or Inf is generated in the testing process, it is
% replaced using the nanapprox function. If there are more than
% NAN_ACCEPT_RATIO of such replacemetes (relative to the size(data,1)), a
% warning message is displayed.
%
% References:
% [Genest and Favre, 2007] Genest, C. and Favre, A. (2007). Everything you always
% wanted to know about copula modeling but were afraid to ask. Hydrol. Eng.,
% 12:347-368.
%
%
% Copyright 2017 Jan Górecki

NAN_ACCEPT_RATIO = 0.05;

switch gofTestName
    case 'K' % Cramér-von Mises statistic based on the Kendall transform 
        Sn = gof2SnK(U, family, theta, emp2copula, NAN_ACCEPT_RATIO);
    case 'R' % Cramér-von Mises statistic based on the Rosenblatt transform 
        Sn = gof2SnR(U, family, theta, NAN_ACCEPT_RATIO); % a precomputed emp2copula wouldn't have sense here 
    case 'E' % Cramér-von Mises statistic based just on the empirical copula
        Sn = gof2SnE(U, family, theta, emp2copula, NAN_ACCEPT_RATIO);
    otherwise 
        error('gofSn: unsupported gof test.');
end
end

%--------------------------------------------------------------------------

function Sn = gof2SnE(data, family, parameter, empCopula, NAN_ACCEPT_RATIO)
%returns Cramer-von Mises statistics based on empirical copula

psiinv = getgenerator(family, parameter, 1);
psi = getgenerator(family, parameter, 0);
yTheo = psi(psiinv(data(:,1))+psiinv(data(:,2)));

% NaN check
[yTheo, nNaNs] = nanapprox(yTheo, data);
if nNaNs/size(data,1) > NAN_ACCEPT_RATIO
    warning('gof2SnE: %d NaNs detected and replaced by their approximations.\n', nNaNs);
end

yEmp = empCopula;

%compute the statistics
Sn = sum((yTheo - yEmp).^2);
end

%--------------------------------------------------------------------------
function Sn = gof2SnK(data, family, parameter, W, NAN_ACCEPT_RATIO)
% Cramér - von Mises statistics for goodness-of-fit test based on the Kendall
% transformation
% for bivariate AC's only
% W - empitical copula of data "W = C(U1,U2)" (precomputed using computeallemp2copulas)

n = size(data,1);

% Kn computation (Kn is the ecdf of W)
Kn = zeros(1,n-1);
for i = 1:(n-1)
    Kn(i) = sum(W <= (i/n));
end
Kn = Kn/n;

% computation of the theoretical distribution W 
% compute Bivariate probability integral transform
% the formulas below are computed using:
% syms t;
% syms theta;
% family = 'A';
% simple(t - getsymbgenerator(family, 1) / diff(getsymbgenerator(family, 1), t))

t = ((1/n):(1/n):1)';
theta = parameter;
switch family
    case 'A'
        K_theta_n = t - (t.*log((t.*theta - theta + 1)./t).*(t.*theta - theta + 1))./(theta - 1);
    case 'C'
        K_theta_n = t - (t.*(t.^theta - 1))./theta;
    case 'F'
        K_theta_n = t + (exp(t.*theta).*log((exp(-t.*theta) - 1)./(exp(-theta) - 1)).*(exp(-t.*theta) - 1))./theta;
    case 'G'
        K_theta_n = t - (t.*log(t))./theta;
    case 'J'
        K_theta_n = t + (log(1 - (1 - t).^theta).*((1 - t).^theta - 1).*(1 - t).^(1 - theta))./theta;
    case '12'
        K_theta_n = t - (t.*(t - 1))./theta;
    case '14'
        K_theta_n = 2.*t - t.*t.^(1./theta);
    case '19'
        K_theta_n = t - (t.^2.*(exp((theta.*(t - 1))./t) - 1))./theta;
    case '20'
        K_theta_n = t - ((3060513257434037.*t.^(theta + 1).*exp(-1./t.^theta))./1125899906842624 - t.^(theta + 1))./theta;
    otherwise
        error('gof2SnK: Unsupported family.');
end
      
% check NaNs
% if a NaN is detected in K_theta, it is replaced by 0, if i < n/2, 
% and by 1, otherwise
nans = isnan(K_theta_n);
nNans = sum(nans);
if (nNans > 0) 
    for i = find(nans)
        if i < n/2
            K_theta_n(i) = 0;
        else
            K_theta_n(i) = 1;
        end
    end
    if (nNans/n > NAN_ACCEPT_RATIO)
        warning('gof2SnK: %d NaNs detected and replaced by their approximations. \n', nNans);
    end
end

% Cramer - von Mises Statistics computation for data
K_theta_n_lin = 0;
K_theta_n_sqr = 0;
for j = 1:n-1
    K_theta_n_lin = K_theta_n_lin + Kn(j)^2 * (K_theta_n(j+1) - K_theta_n(j));
    K_theta_n_sqr = K_theta_n_sqr + Kn(j) * (K_theta_n(j+1)^2 - K_theta_n(j)^2);
end 
Sn = n/3 + n * K_theta_n_lin - n * K_theta_n_sqr;
end

%--------------------------------------------------------------------------

function Sn = gof2SnR(data, family, parameter, NAN_ACCEPT_RATIO)
% return the S_n^{(C)} (based on the Rosenblatt's transformation)
%
% NOTE: if no real value is returned after the transformation, the original
% parameter is divided by 10 and the transformation is re-computed until
% some real value is obtained.

u1 = data(:,1);
u2 = data(:,2);
theta3 = parameter;

isAtLeastOneNotNaN = false;
while ~(isAtLeastOneNotNaN)
    
    % \delta C / \delta u1 (u1, u2, theta3)
    % the formulas are obtained through:
    % 0) syms u1; syms u2; syms theta3;
    % 1) AC = HACopula({{'20', 1.5}, 1, 2});
    % 2) simple(diff(getcdf(AC), u1))
    
    switch family
        case 'A'
            e2 = (u2 + theta3.*u2.*(u2 - 1))./(theta3.*u2 - theta3 + u1.*(theta3 - theta3.*u2) + 1).^2;
        case 'C'
            e2 = 1./(u1.^(theta3 + 1).*(1./u1.^theta3 + 1./u2.^theta3 - 1).^(1./theta3 + 1));
        case 'F'
            e2 = (exp(theta3.*(u2 + 1)) - exp(theta3))./(exp(theta3.*(u1 + 1)) + exp(theta3.*(u2 + 1)) - exp(theta3.*(u1 + u2)) - exp(theta3));
        case 'G'
            e2 = (exp(-((-log(u1)).^theta3 + (-log(u2)).^theta3).^(1./theta3)).*(-log(u1)).^(theta3 - 1).*((-log(u1)).^theta3 + (-log(u2)).^theta3).^(1./theta3 - 1))./u1;
        case 'J'
            e2 = -((1 - u2).^theta3 - 1).*(1 - u1).^(theta3 - 1).*((1 - u1).^theta3 + (1 - u2).^theta3 - (1 - u1).^theta3.*(1 - u2).^theta3).^(1./theta3 - 1);
        case '12'
            e2 = ((1./u1 - 1).^(theta3 - 1).*((1./u1 - 1).^theta3 + (1./u2 - 1).^theta3).^(1./theta3 - 1))./(u1.^2.*(((1./u1 - 1).^theta3 + (1./u2 - 1).^theta3).^(1./theta3) + 1).^2);
        case '14'
            e2 = ((1./u1.^(1./theta3) - 1).^(theta3 - 1).*((1./u1.^(1./theta3) - 1).^theta3 + (1./u2.^(1./theta3) - 1).^theta3).^(1./theta3 - 1))./(u1.^(1./theta3 + 1).*(((1./u1.^(1./theta3) - 1).^theta3 + (1./u2.^(1./theta3) - 1).^theta3).^(1./theta3) + 1).^(theta3 + 1));
        case '19'
            e2 = (theta3.^2.*exp(theta3./u1))./(u1.^2.*log(exp(theta3./u1) + exp(theta3./u2) - exp(theta3)).^2.*(exp(theta3./u1) + exp(theta3./u2) - exp(theta3)));
        case '20'
            e2 = exp(1./u1.^theta3)./(u1.^(theta3 + 1).*log(exp(1./u1.^theta3) + exp(1./u2.^theta3) - exp(1)).^(1./theta3 + 1).*(exp(1./u1.^theta3) + exp(1./u2.^theta3) - exp(1)));
        otherwise
            error('gof2SnR: Unsupported family');
    end
    
    % provide at least one real number in e2
    if sum(~isnan(e2)) == 0
        % there is no real number in e2

        % decrease the parameter 
        theta3 = theta3/10;
        
        % check the range
        famRange = tau2theta(family, getfamilytaurange(family));
        if theta3 < famRange(1)
            error('gof2SnR: Unable to find a parameter value such that the Rosenblatt''s transformation would return at least one real value.');
        else
            warning(['gof2SnR: the Rosenblatt''s transformation has returned for the original parameter ' num2str(parameter) ' only NaN for all data. Trying the parameter' num2str(theta3) ' instead.']);
        end
    else
        isAtLeastOneNotNaN = true;
    end
        
end

% NaN check
[e2, nNaNs] = nanapprox(e2, data);
if nNaNs/size(data,1) > NAN_ACCEPT_RATIO
    warning('gof2SnR: %d NaNs detected and replaced by their approximations.\n', nNaNs);
end

%Sn(C) statistics according to Genest
e1 = u1;
C_pi = e1 .* e2;

% D computation
D = computeallemp2copulas([e1 e2]);
D = D{1,2};

Sn = sum((D - C_pi).^2);
end




