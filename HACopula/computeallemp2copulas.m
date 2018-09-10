function emp2copulas = computeallemp2copulas(U)
%COMPUTEALLEMP2COPULAS - compute all bivariate empirical copulas.
%
% emp2copulas = computeallemp2copulas(U) returns a d x d cell array emp2copulas,
% where emp2copulas{i,j} contains a vector of size n (n =
% size(U,1)) describing the bivariate empirical copula computed from 
% i-th and j-th columns of U.
%
% NOTE:
% As emp2copulas{i,j} == emp2copulas{j,i} for all i,j from {1, ..., d} (d =
% size(U,2)), all cells with i > j are empty.
%
% References:
% [Genest and Favre, 2007] Genest, C. and Favre, A. (2007). Everything you always
% wanted to know about copula modeling but were afraid to ask. Hydrol. Eng.,
% 12:347-368.
%
%
% Copyright 2018 Jan Gorecki

[n, d] = size(U);

emp2copulas = cell(d);

for i = 1:d
    for j = i + 1:d
        emp2copulas{i,j} = zeros(n,1);
        for ii=1:n
            emp2copulas{i,j}(ii) = sum((U(:,i) <= U(ii,i)) .* (U(:,j) <= U(ii,j)));
        end
        emp2copulas{i,j} = emp2copulas{i,j} / n;
    end
end

end
