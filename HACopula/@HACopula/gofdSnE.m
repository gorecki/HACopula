function Sn = gofdSnE(obj, U)
%GOFDSNE - returns the value of the Cramer-von Mises statistics (denoted
% S_n in [Genest et al., 2009]) based on the empirical copula for
% the HAC obj and the observations U.
%
%
% Copyright 2018 Jan Gorecki

% perform basic data checks
if ~iscopuladata(U)
    warning('HACopula:notCopulaData', 'HACopula::gofdSnE: At least one column of U was rejected to be standard uniform according to KS test at the 5 %% significance level, i.e., U might not be a sample from a copula.');
end

%compute the "empirical copula" for the data
[n, d] = size(U);

Cn = zeros(n,1);
for i=1:n
    mult = (U(:,1) <= U(i,1));
    for j=2:d
        mult = mult .* (U(:,j) <= U(i,j));
    end
    Cn(i) = sum(mult);
end
Cn = Cn / n;

% slow symbolic evaluation
%get theoretical cdf of the copula for estimated parameter
%             cdfCtheta = getcdf(obj);

%compute statistics
%             forksArray = getforksarray(obj);
%             Usymb = sym('u%d', [1 d]);
%             theta = sym('theta%d', [1 2*d-1]);
%Ctheta = evalsymbformula(cdfCtheta, [Usymb, theta([forksArray(:).TauOrdering])], U, [forksArray(:).Parameter]);

% fast non-symbolic evaluation
Ctheta = cdf(obj, U);

% NaN and Inf check
[Ctheta, nNaNs] = nanapprox(Ctheta, U);
if nNaNs > 0
    warning('HACopula:NaN_detected', ['HACopula:gofdSnE:: ' num2str(nNaNs) ' NaNs detected and replaced by their approximations.']);
end

Sn = sum((Ctheta - Cn).^2);
end