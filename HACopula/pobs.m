function U = pobs(X, method)
%POBS - Compute uniform pseudo-observations for a random sample
% Inputs:
% X         - a random sample
% method    - 'rank' (default, if *method* is not provided) or 'kernel'.
%             For method == 'rank', the pseudo-observations are computed
%             according to the approach described in [Genest and Favre,
%             2007], where n+1 is used for division to keep the
%             pseudo-observations lower than 1.
%             For method == 'kernel', the pseudo-observations are computed
%             using Kernel smoothing function estimate for univariate data
%             d-times. This method is useful in a case when one wants to
%             trasform from the pseudo-observations (U) back to the original
%             observations (X), e.g., in Bayesian classfication, see
%             [Górecki et al., 2016a].
%             
%
% [Genest and Favre, 2007] Genest, C. and Favre, A. (2007). Everything you always
%     wanted to know about copula modeling but were afraid to ask. Hydrol. Eng.,
%     12:347-368.
% [Górecki et al., 2016a] Górecki, J., Hofert, M., and Holeòa, M. (2016). An 
%     approach to structure determination and estimation of hierarchical
%     Archimedean copulas and its application to bayesian classication.
%     Journal of Intelligent Information Systems, pages 21-59.
%
% Copyright 2017 Jan Górecki

if nargin == 1
    method = 'rank';
elseif nargin ~= 2
    error('pobs: Invalid number of input arguments.');
end

[n, d] = size(X);
U = zeros(n, d);

for i=1:d
    switch method
        case 'rank'
            U(:,i) = getmaxrank(X(:,i)) / (n + 1);
        case 'kernel'
            U(:,i) = ksdensity(X(:,i), X(:,i), 'function', 'cdf');
        otherwise
            error(['pobs: ' method ' is unknown method. Choose one from {''rank'', ''kernel''}']);
    end
end

end

function maxRank = getmaxrank(X)
n = size(X, 1);
maxRank = zeros(n, 1);
[sorted, I] =  sort(X);
r = n;
prev = sorted(n);

% proceed in the reverse order
for i=n:-1:1
    x = sorted(i);
    if x == prev
        maxRank(I(i)) = r;
    else
        prev = x;
        r = i;
        maxRank(I(i)) = r;
    end    
end

end