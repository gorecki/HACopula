function pvalue = computepvalue(obj, U, estimator, N)
% COMPUTEPVALUE - bootstrap computation of p-value
%
% Purpose:
% Given a HAC estimate, pseoudo-observations and an estimator,
% it performs the bootstrap method proposed in [Genest et al.,
% 2009] in order to get a p-value estimate.
%
% Inputs:
% obj       - A HAC estimate obtained by the input estimator
%             for the input data U.
% U         - Pseudo-observations.
% estimator - A HAC estimator with only one input argument
%             corresponding to pseudoobservations, e.g., an
%             anonymous funtion @(U)f(U).
% N         - A desired number of bootstrap replications.
%
% Outputs:
% pvalue    - A value from [0, 1], which can be interpreted as
%             the probability (under valid H0, i.e., the
%             underlying copula is obj) that, using a given HAC
%             estimator, the goodness-of-fit, would be worse
%             than the actual observed results, i.e.,
%             goodness-of-fit of the
%             estimate obj (obtained by the input estimator)
%             for the data U.
%
% NOTE:
% The goodness-of-fit statistic used in the computation is the
% statistic denoted S_n in [Genest et al., 2009].
%
%
% Copyright 2017 Jan Górecki

n = size(U,1);
estSn = gofdSnE(obj, U);
[~,theirSn] = bootstrapSn(obj, estimator, N, n);
pvalue = sum(theirSn > estSn)/size(theirSn,2);
end