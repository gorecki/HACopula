function [bootstrapSamples,theirSn] = bootstrapSn(obj, estimator, bootstrapSize, sampleSize)
% BOOTSTRAPSN - a bootstrap for a copula estimate
%
% Inputs:
% obj           - A HACopula object.
% estimator     - A HAC estimator with one input argument
%                 corresponding to pseudoobservations, e.g., an
%                 anonymous funtion @(U)f(U).
% bootstrapSize - A bootstrap size (default 1000).
% sampleSize    - A sample size (default 100).
%
% For a given HAC, performs the bootstrap method proposed in
% [Genest et al., 2009], returns the bootstrap samples and
% their values of the rank-based Cramer-von Mises statistics
% (denoted Sn in the paper).
%
%
% Copyright 2018 Jan Gorecki

narginchk(2,4);

if nargin < 3
    bootstrapSize = 1000;
    sampleSize = 100;
elseif nargin < 4
    sampleSize = 100;
end

% start bootstrapping process

bootstrapSamples = cell(1, bootstrapSize);
theirSn = zeros(1, bootstrapSize);

for sample = 1:bootstrapSize
    %                 if mod(sample,10) == 1
    %                     disp(sample);
    %                 end
    % sample
    UKnown = rnd(obj, sampleSize);
    % turn to pseudo-observations
    U = pobs(UKnown);
    bootstrapSamples{sample} = U;
    % compute GoF for the estimate on the sample
    theirSn(sample) = gofdSnE(estimator(U), U);
end
end
