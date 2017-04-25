function out = evalsurv(obj, U)
%EVALSURV - evaluates the survival copula of obj at U \in [0, 1]^d
%
% NOTE:
% The computation involves evaluations of the HAC obj in all 2^d
% corners of the d-dimensional hypercube (U, 1).
%
%
% Copyright 2017 Jan Górecki

out = prob(obj, U, 1*ones(1, getdimension(obj)));
end