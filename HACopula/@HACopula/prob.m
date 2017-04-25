function out = prob(obj, l, u)
% PROP - probability of falling into a hypercube from [0,1]^d
%
% Compute the probability of a d-dimensional random vector U
%   distributed according to a given copula *obj* to fall in a
%   hypercube (l, u], where *l* and *u* denote the lower and upper
%   corners of the hypercube, respectively.
%
% NOTE:
% The computation involves evaluations of the HAC obj in all 2^d
% corners of the hypercube.
%
%
% Copyright 2017 Jan Górecki

d = getdimension(obj);
% get the corners of the d-dimensional unit hypercube.
[corner{1:d}] = ndgrid(logical([0 1]));
corner = cat(d+1,corner{d:-1:1});
corner = double(reshape(corner,[],d));

% get the corners of [l u]^d
lu = zeros(size(corner));
for i = 1:d
    lu(:, i) = corner(:, i) * (u(i)-l(i)) + l(i);
end

% evaluate the copula at these corners
cornProb = evaluate(obj, lu);

% compute the C-volume of [l u]^d
out = 0;
for i = 1:size(corner, 1)
    out = out + (-1)^sum(1-corner(i,:)) * cornProb(i);
end
end