function pValue = computepvalue(obj, U, estimator, N)
% COMPUTEPVALUE - bootstrap computation of p-value
%
% Purpose:
% Given a HAC estimate, pseoudo-observations and an estimator,
% it performs the bootstrap method proposed in [Genest et al.,
% 2009] in order to get a p-value estimate.
%
% IMPORTANT:
% Even if the implementation of the parametric bootstrap method exactly
% follows its theoretical description proposed in [Genest et al., 2009],
% the requirements on the properties of the involved rank-based estimator
% addressed in [Genest and Remillard, 2008] in order to guarantee that the
% parametric bootstrap method yields valid p-values for S_n haven been
% theoretically proven yet.  More precisely, according to [Genest and Remillard, 2008],
% the mentioned requirements are known to be satisfied only for 2-ACs,
% i.e., just for a sub-class of all HACs supported by HACopula. Any
% results obtained using this method thus should be interpreted
% with this in mind.
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
% References:
% [Genest and Remillard, 2008] Genest, C., Remillard, B. (2008). Validity
%   of the parametric bootstrap for goodness-of-fit: testing in
%   semiparametric models. In Annales de lÿInstitut Henri Poincare:
%   Probabilites et Statistiques, volume 44, pp. 1096ÿ1127
%
% [Genest et al., 2009] Genest, C., Remillard, B., and Beaudoin, D.
%    (2009). Goodness-of-fit tests for copulas: A review and a power
%    study. Insurance: Mathematics and Economics, 44(2):199-213.
%
%
%
% Copyright 2018 Jan Gorecki

n = size(U,1);
estSn = gofdSnE(obj, U);
[~,theirSn] = bootstrapSn(obj, estimator, N, n);
pValue = sum(theirSn > estSn)/size(theirSn,2);
end