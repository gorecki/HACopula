function out = cdf(obj, U)
%CDF - evaluates CDF of the HAC obj at a certain point (U)
%
% Assuming that obj represents a copula function C, the output
% is C(U), where C is applied row-wise on the rows of U.
%
%
% Copyright 2018 Jan Gorecki

n = size(U, 1);
nChildren = size(obj.Child,2);
subU = zeros(n, nChildren);
for i = 1:nChildren
    if isa(obj.Child{i}, 'HACopula')
        % go to recursion
        subU(:,i) = cdf(obj.Child{i}, U);
    else % it is a leaf
        % return the column corresponding to the leaf
        subU(:,i) = U(:,obj.Child{i});
    end
end
psi = getgenerator(obj.Family, obj.Parameter, 0);
psiInv = getgenerator(obj.Family, obj.Parameter, 1);
out = psi(sum(psiInv(subU), 2));
end
