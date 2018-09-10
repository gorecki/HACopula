function UTrans = computediagonaltransform(obj, U)
%COMPUTEDIAGONALTRANSFORM - transformes U(:, obj.Leaves) according to obj
% into a single vector of observations through the diagonal
% transformation (see [Gorecki et al., 2014]).
%
% TODO: use this method in collapse for MLE in the diagonal
% estimation
%
%
% Copyright 2018 Jan Gorecki

% prepare necessary functions
psi = getgenerator(obj.Family, obj.Parameter, 0);
psiInv = getgenerator(obj.Family, obj.Parameter, 1);
nChild = length(obj.Child);
delta = @(t) psi(nChild.* psiInv(t));  % the transformation

% prepare necessary data
UChild = zeros(size(U,1), nChild);
for i = 1:nChild
    if isa(obj.Child{i}, 'HACopula')
        % go into recursion
        UChild(:, i) = computediagonaltransform(obj.Child{i}, U);
    else
        % use the observations from U
        UChild(:, i) = U(:, obj.Child{i});
    end
end

% compute the transformation
UTrans = delta(max(UChild,[], 2));

% do a NaN check
[UTrans, nNaNs] = nanapprox(UTrans, UChild);
if nNaNs > 0
    warning('HACopula:NaN_detected', ['HACopula::computediagonaltransform:: There was ' num2str(nNaNs) ' NaNs detected in the diagonal transformation and replaced by their approximations.']);
end
end