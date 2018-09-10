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
% NAN_ACCEPT_RATIO of such replacemetes (relative to size(data,1)), a
% warning message is displayed.
%
% References:
% [Genest and Favre, 2007] Genest, C. and Favre, A. (2007). Everything you always
% wanted to know about copula modeling but were afraid to ask. Hydrol. Eng.,
% 12:347-368.
%
%
% Copyright 2018 Jan Gorecki

NAN_ACCEPT_RATIO = 0.05;

switch gofTestName
    case 'K' % Cramer-von Mises statistic based on the Kendall transform 
        Sn = gof2SnK(U, family, theta, emp2copula, NAN_ACCEPT_RATIO);
    case 'R' % Cramer-von Mises statistic based on the Rosenblatt transform 
        Sn = gof2SnR(U, family, theta, NAN_ACCEPT_RATIO); % a precomputed emp2copula wouldn't have sense here 
    case 'E' % Cramer-von Mises statistic based just on the empirical copula
        Sn = gof2SnE(U, family, theta, emp2copula, NAN_ACCEPT_RATIO);
    otherwise 
        error('HACopula:BadInputs', 'gofSn: unsupported gof test.');
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
    warning('HACopula:NaN_detected', 'gof2SnE: %d NaNs detected and replaced by their approximations.\n', nNaNs);
end

yEmp = empCopula;

%compute the statistics
Sn = sum((yTheo - yEmp).^2);
end

%--------------------------------------------------------------------------
function Sn = gof2SnK(data, family, parameter, W, NAN_ACCEPT_RATIO)
% Cramer - von Mises statistics for goodness-of-fit test based on the Kendall
% transformation
% for bivariate AC's only
% W - empirical copula of data "W = C(U1,U2)" (precomputed using computeallemp2copulas)

n = size(data,1);

% Kn computation (Kn is the ecdf of W)
Kn = zeros(1,n-1);
for i = 1:(n-1)
    Kn(i) = sum(W <= (i/n));
end
Kn = Kn/n;

% computation of the theoretical distribution W 
K_theta_n = ACkendtrans(family, parameter, ((1/n):(1/n):1)');
      
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
        warning('HACopula:NaN_detected', ['gof2SnK: ' num2str(nNans) ' NaNs detected and replaced by their approximations.']);
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
theta = parameter;

isAtLeastOneNotNaN = false;
while ~(isAtLeastOneNotNaN)
    
  
    e2 = ACcondcdf(family, theta, u1, u2);
    
    % provide at least one real number in e2
    if sum(~isnan(e2)) == 0
        % there is no real number in e2

        % decrease the parameter 
        theta = theta/10;
        
        % check the range
        famRange = tau2theta(family, getfamilytaurange(family));
        if theta < famRange(1)
            error('HACopula:gof2Sn', 'gof2SnR: Unable to find a parameter value such that the Rosenblatt''s transformation would return at least one real value.');
        else
            warning('HACopula:gof2Sn', ['gof2SnR: the Rosenblatt''s transformation has returned for the original parameter ' num2str(parameter) ' only NaN for all data. Trying the parameter' num2str(theta) ' instead.']);
        end
    else
        isAtLeastOneNotNaN = true;
    end
        
end

% NaN check
[e2, nNaNs] = nanapprox(e2, data);
if nNaNs/size(data,1) > NAN_ACCEPT_RATIO
    warning('HACopula:NaN_detected', ['gof2SnR: ' num2str(nNaNs) ' NaNs detected and replaced by their approximations.']);
end

%Sn(C) statistics according to Genest
e1 = u1;
C_pi = e1 .* e2;

% D computation
D = computeallemp2copulas([e1 e2]);
D = D{1,2};

Sn = sum((D - C_pi).^2);
end




