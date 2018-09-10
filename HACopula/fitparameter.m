function theta = fitparameter(descLeaves, family, hacEstimatorRatio, thetaEstimatorPairwise, ...
                               thetaEstimatorDiagonal, K, U, g1, iK, jK, maxSimilarity, r, attitude)
%FITPARAMETER - Fitting the parameter of an Archimedean generator
%
% theta = fitparameter(descLeaves, family, hacEstimatorRatio,
% thetaEstimatorPairwise, thetaEstimatorDiagonal, K, U, g1, iK, jK,
% maxSimilarity, r) returns an estimate of the generator such that all its
% children have thier descendant leaves stored in the input descLeaves.
%
% Inputs:
% descLeaves            - The descendant leaves of the generator for which
%                         the parameter should be estimated, e.g., if {[1 2],
%                         [3 4]} is supplied, then U13, U14, U23 and U24 are considered.
% family                - A family assumed for the generator.
% hacEstimatorRatio     - A value from [0, 1].
% thetaEstimatorPairwise- An estimator for the parameter of ACs (see below for details).
% thetaEstimatorDiagonal- An estimator for the parameter of ACs (see below for details).
% K                     - The matrix of Kendall's tau coefficients for U.
% U                     - Pseudo-observations.
% g1                    - An [0, 1] aggregation function (g in Algorithm 3 from [Gorecki et al., 2016a]).
% iK, jK                - The indices of the columns in U from which the parameter
%                         will be estimated (related to the diagonal HAC
%                         estimation).
% maxSimilarity         - A value from [0, 1] (tau-based-estimate of the parameter).
% r                     - An admissible range for the parameter (used for MLE estimation).
% attitude              - 'optimistic' or 'pessimistic'.
%
% Outputs:
% theta         - An estimate of the parameter.
%
%
%
% NOTE:
% If no admissible estimate is found, theta = NaN,
%                 otherwise, theta is a real number
% If the HAC estimator is chosen 'pairwise' or combined (hacEstimatorRatio < 1), one out of the 
% 4 different AC estimators can be chosen:
%
% invtau    - it first aggregates corresponding taus through g1 and then inverts it
%             through tau^{-1}_{family}
% invtau2   - it first inverts all corresponding taus through
%             tau^{-1}_{family} and then aggregates it through g1
% mle       - it aggregates all corresponding bivariate margins, e.g., if
%             (U1, U2) and (U1, U3) corresponds to the generator, these margins are
%             aggregated to [U1 U2; U1 U3] and used in the MLE
% mle2      - it first estimates all thetas for corresponding bivariate
%             margins using the MLE and then aggregates these thetas using g1
%
% If the HAC estimator is chosen 'diagonal' or combined (hacEstimatorRatio > 0), then there 
% are only 2 different AC estimators, because invtau is the same as invtau2 
% and mle is the same as mle2 due to the diagonal HAC estimation approach:
%
% invtau (= invtau2)    - inverts the corresponding tau estimate through
%                         tau^{-1}_{family}
% mle (= mle2)          - uses the MLE for the corresponding bivariate
%                         margin
%
% If the HAC estimator is chosen to be a mix of pairwise and diagonal (0 <
% hacEstimatorRatio < 1), the AC parameter estimate is a weighted (by
% hacEstimatorRatio) average of the pairwise and diagonal estimates.
%
% References:
% [Gorecki et al., 2014] Gorecki, J., Hofert, M., and Holena, M. (2014). On
%     the consistency of an estimator for hierarchical Archimedean copulas.
%     In 32nd International Conference on Mathematical Methods in
%     Economics, pages 239-244.
% [Gorecki et al., 2016a] Gorecki, J., Hofert, M., and Holena, M. (2016). An 
%     approach to structure determination and estimation of hierarchical
%     Archimedean copulas and its application to bayesian classication.
%     Journal of Intelligent Information Systems, pages 21-59.
% [Gorecki et al., 2017] On Structure, Family and Parameter Estimation
%     of Hierarchical Archimedean copulas. Journal of Statistical Computation 
%     and Simulation, 87(17), 3261-3324
%
%
% Copyright 2018 Jan Gorecki

%--------------------------------------------------------------------------


% assume a NaN for the estimates 
thetaActPairwise = NaN;
thetaActDiagonal = NaN;

nChildren = size(descLeaves,2);
if hacEstimatorRatio < 1 % pairwise or combined
    switch thetaEstimatorPairwise
        case 'invtau' % first aggregate then invert tau to theta
            vecSubK = [];
            % aggregate all margins
            for i = 1:nChildren
                for j = i+1:nChildren
                    subK = K(descLeaves{i}, descLeaves{j});
                    vecSubK = [vecSubK reshape(subK, [1, length(descLeaves{i})*length(descLeaves{j})])];
                end
            end
            tauEst = g1(vecSubK);
            
            if strcmp(attitude, 'optimistic')
                % set the parameter to the lower (or upper) bound of the
                % parameter range in order to get at least some estimate (even
                % if low-fitting)
                famRange = getfamilytaurange(family);
                if tauEst < famRange(1)
                    tauEst = famRange(1);
                end
                if tauEst > famRange(2)
                    tauEst = famRange(2);
                end
            end
            
            % now assume both attitudes
            if istauinvertible(family, tauEst)
                thetaActPairwise = tau2theta(family, tauEst);
            else
                thetaActPairwise = NaN;
            end
        case 'invtau2' % first invert taus and then aggreagate
            vecSubK = [];
            % aggregate all margins
            for i = 1:nChildren
                for j = i+1:nChildren
                    subK = K(descLeaves{i}, descLeaves{j});
                    vecSubK = [vecSubK reshape(subK, [1, length(descLeaves{i})*length(descLeaves{j})])];
                end
            end
            
            if strcmp(attitude, 'optimistic')
                famRange = getfamilytaurange(family);
                % trim all particular estimates to invertable ranges
                for m = 1:size(vecSubK, 2)
                    % set the parameter to the lower (or upper) bound of the
                    % parameter range in order to get at least some estimate (even
                    % if low-fitting)
                    if vecSubK(m) < famRange(1)
                        vecSubK(m) = famRange(1);
                    end
                    if vecSubK(m) > famRange(2)
                        vecSubK(m) = famRange(2);
                    end
                end
            end
            
            % now assume both attitudes
            if istauinvertible(family, vecSubK)
                thetaActPairwise = g1(tau2theta(family, vecSubK));
            else
                % try tauinv1 if tauinv2 does not return a number
                tauEst = g1(vecSubK);

                if istauinvertible(family, tauEst)
                    thetaActPairwise = tau2theta(family, tauEst);
                else
                    thetaActPairwise = NaN;
                end
            end
            
        case 'mle'  
            % aggregate all bimarginal data corresponding to the node, e.g., [U1 U2; U1 U3], and use it for the MLE
            % NOTE: faster than mle2
            
            biU = [];
            % aggregate all bimarginal data
            for i = 1:nChildren
                for j = i+1:nChildren
                    for ii = 1:length(descLeaves{i})
                        for jj = 1:length(descLeaves{j})
                            biU = [biU ; U(:,[descLeaves{i}(ii) descLeaves{j}(jj)])];
                        end
                    end
                end
            end
            thetaActPairwise = mlefor2AC(biU, family, maxSimilarity, r(1), r(2));
            if ~((thetaActPairwise > -Inf) && (thetaActPairwise < Inf))
                error('HACopula:NaN_detected', 'fitparameter: MLE has returned a NaN.');
            end
            
        case 'mle2'
            vecMlePair = [];
            % aggregate all margins
            for i = 1:nChildren
                for j = i+1:nChildren
                    mlePair = zeros(length(descLeaves{i}), length(descLeaves{j}));
                    for ii = 1:length(descLeaves{i})
                        for jj = 1:length(descLeaves{j})
                            mlePair(ii, jj) = mlefor2AC(U(:,[descLeaves{i}(ii) descLeaves{j}(jj)]), family, maxSimilarity, r(1), r(2));
                            if ~((mlePair(ii, jj) > -Inf) && (mlePair(ii, jj) < Inf))
                                error('HACopula:NaN_detected', 'fitparameter: MLE has returned a NaN.');
                            end
                        end
                    end
                    
                    % compute the resulting theta as the g1-aggregated value of the
                    % pairwise theta estimatimates
                    vecMlePair = [vecMlePair reshape(mlePair, [1, length(descLeaves{i})*length(descLeaves{j})])];
                end
            end
            thetaActPairwise = g1(vecMlePair);
            
        otherwise
            error('HACopula:BadInputs', 'fitparameter: Unsupported pairwise estimator. Should be ''invtau'', ''invtau2'', ''mle'' or ''mle2''.')
    end
end

% NOTE: diagonal estimation is not used for collapsing (only for
% HAC estimation)
if hacEstimatorRatio > 0 % diagonal or combined
    switch thetaEstimatorDiagonal
        case {'invtau', 'invtau2'} % in the diagonal case, these two estimators are equal
            
            tauEst = K(iK, jK);

            if strcmp(attitude, 'optimistic')
                % set the parameter to the lower (or upper) bound of the
                % parameter range in order to get at least some estimate (even
                % if low-fitting)
                famRange = getfamilytaurange(family);
                if tauEst < famRange(1)
                    tauEst = famRange(1);
                end
                if tauEst > famRange(2)
                    tauEst = famRange(2);
                end
            end
            
            % now assume both attitudes
            if istauinvertible(family, tauEst)
                thetaActDiagonal = tau2theta(family, tauEst);
            else
                thetaActDiagonal = NaN;
            end
            
        case {'mle', 'mle2'}
            thetaActDiagonal = mlefor2AC(U(:,[iK jK]), family, maxSimilarity, r(1), r(2)); 
            if ~((thetaActDiagonal > -Inf) && (thetaActDiagonal < Inf))
                error('HACopula:NaN_detected', 'fitparameter: MLE has returned a NaN.');
            end
        otherwise
            error('HACopula:BadInputs', 'fitparameter: Unsupported diagonal estimator. Should be ''invtau'', ''invtau2'', ''mle'' or ''mle2''.')

    end
end

% is the theta a real number?, i.e., check if a NaN resulted in the
% estimation process
isThetaActPairwiseReal = (thetaActPairwise > -Inf) && (thetaActPairwise < Inf);
isThetaActDiagonalReal = (thetaActDiagonal > -Inf) && (thetaActDiagonal < Inf);

% this is the optimistic version
if (~isThetaActPairwiseReal) && (~isThetaActDiagonalReal)
    theta = NaN;
elseif isThetaActPairwiseReal && (~isThetaActDiagonalReal)
    theta = thetaActPairwise;
elseif isThetaActDiagonalReal && (~isThetaActPairwiseReal)
    theta = thetaActDiagonal;
else
    theta = hacEstimatorRatio * thetaActPairwise + (1-hacEstimatorRatio) * thetaActDiagonal; % weighted average
end

end
