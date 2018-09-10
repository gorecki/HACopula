function [HACObject, fitLog] = HACopulafit(U, families, varargin)
% HACOPULAFIT - Estimate a hierarchical Archimedean copula (HAC)
%
% Purpose:
% The function returns an estimate of the underlying copula assuming it is
% a HAC. Given a matrix of (pseudo-)observations U from [0, 1] and a set of
% Archimedean families, it returns an instance of the HACopula class, i.e.,
% a d-HAC model, which has all the generators from the given set of
% families. If #families > 1, then it might return either a heterogenous HAC or
% a homogeneous HAC. If #families == 1, it always returns a homogenous HAC. 
%
% Usage:
% HACObject = HACopulafit(U, families)
% - returns a HAC estimate
% [HACObject, fitLog] = HACopulafit(U, families)
% - returns a HAC estimate and a log of the estimation process
% 
%
% Inputs (obligatory):
% U                 - An (n * d)-matrix of (pseudo-)observations from [0, 1]
%                     Note that each element of U is checked to be from [0,
%                     1] and also each column is checked to be uniformly
%                     distributed using two-sample Kolmogorov-Smirnov test
%                     against a sample with a perfect standard uniform
%                     distribution for all one-dimensional margins. These
%                     checks could be switched off using the parameter
%                     'CheckData'.
% families          - A cell containing labels of Archimedean families of 
%                     generators. The possible labels are: 'A', 'C', 'F',
%                     'G', 'J', '12', '14', '19', '20', '?', where the
%                     latter corresponds to an arbitary 'unknown' family,
%                     which serves just for purposes of structure
%                     estimation and collapsing under no assumptions on the
%                     underlying families, and thus cannot be combined with
%                     other families. For heterogenous estimation, any
%                     subset of either {'A', 'C', '19', '20'} or {'C',
%                     '12', '14', '19', '20'} is possible. Note that for
%                     'F', 'G' and 'J', only homogeneous estimation is
%                     possible. Example: estimate = HACopulafit(U, {'C',
%                     '19'})
%
% Inputs (optional):
% HACEstimator      - If 'HACEstimator' == 'pairwise' (or  == 0), the HAC is
%                     estimated according to the approach represented by
%                     Algorithms 2 and 3 in [Gorecki et al., 2017]. If
%                     'HACEstimator' == 'diagonal' (or  == 1), the HAC is
%                     estimated according to the approach from [Gorecki et
%                     al., 2016b] represented by Algorithm 4. If
%                     'HACEstimator' is in (0, 1), both of the approaches
%                     are used at the same time and their results are
%                     weighted according to the value of the parameter,
%                     e.g., theta = (1-HACEstimator) * thetaPairwise +
%                     HACEstimator * thetaDiagonal. Note that the diagonal
%                     approach is independent both on the aggregation
%                     functions g_1 and g_2.
% ThetaEstimator    - An estimator for the parameter of an AC generator in the
%                     estimated HAC for the *pairwise* HAC estimation. The
%                     options are 'invtau', 'invtau2', 'mle' and 'mle2'.
% ThetaEstimator2   - An estimator for the parameter of an AC generator in the
%                     estimated HAC for the *diagonal* HAC estimation. The
%                     options are 'invtau', 'invtau2', 'mle' and 'mle2',
%                     however, 'invtau' = 'invtau2' and 'mle' = 'mle2'.
%                     Note that if ThetaEstimator2 is not specified, it is
%                     set to the same value as ThetaEstimator. 
% g_1               - A [0, 1]-aggregation function for aggregation of theta
%                     estimates from bivariate margins. The options are
%                     'average', 'min' or 'max'.
% g_2               - A [0, Inf)-aggregation function for aggregation of gof
%                     statistics from bivariate margins. 'g_2' could be any
%                     [0, Inf)-aggregation function encoded as an anonymous
%                     function, e.g., @(t)max(t). Note that if #families ==
%                     1, then g_2 does not influence the result. 
% GOF               - If #families > 1, then a GOF statistic is used to
%                     select an appropriate family for an estimated
%                     generator. The options are 'E' refering to S_n^{(E)},
%                     'K' refering to S_n^{(K)} (Kendall's transform) and 'R'
%                     refering to S_n^{(C)} (Rosenblatt's transform) from
%                     [Gorecki et al., 2017]. Note that the chosen statistic
%                     is first computed for all bivariate margins and then
%                     aggregated with g_2.
% KendallMatrix     - The Kendall correlation matrix of U. Use (could be
%                     pre-computed by kendallTauMatrix(U)), if more
%                     estimation processes are done with the same U.
% Emp2copulas       - The bivariate empirical copulas precomputed from U.
%                     Use (could be pre-computed by computeallemp2copulas(U)),
%                     if doing the heterogenous estimation and more
%                     estimation processes are done with the same U.
% Attitude          - If 'Attitude' == 'pessimistic', it is 
%                     possible, if no appropriate estimate could be found,
%                     e.g., for the family labeled A, if tau >= 1/3 and
%                     'ThetaEstimator' = 'invtau', that no HAC is returned,
%                     precisely, [] is returned. Otherwise, i.e, if
%                     'Attitude' == 'optimistic', some estimate is always
%                     returned (parameters are truncated in order to get a
%                     proper copula).
% CheckData         - If 'CheckData' == 'on', each element of U is checked 
%                     to be from [0, 1] and also each column is checked to
%                     be uniformly distributed using the two-sample
%                     Kolmogorov-Smirnov test against a sample with a
%                     perfect standard uniform distribution for all
%                     one-dimensional margins. If 'CheckData' == 'off', no
%                     checks on U are done.
% PreCollapsedHAC   - A HAC (an instance of the HACopula class) that has
%                     been already collapsed to a non-binary one. Providing
%                     such a HAC, no structure estimation is performed and
%                     the structure of PreCollapsedHAC is used instead.
%                     Based on this structure, HACopulafit estimates only
%                     the parameters and the families.
% PreCollapse       - If == true (possible only for HACEstimator == 0, i.e., pairwise HAC estimator), 
%                     HACopulafit first estimates a binary HAC structure with no
%                     assumptions on the underlying families (i.e.,
%                     internally using the family '?'), then collapses to a
%                     non-binary one and finally supply the resulting
%                     structure to HACopulafit instead of PreCollapsedHAC.
%                     If == false, no pre-collapsing is performed. Note
%                     that in the heterogeneous case, pre-collapsing shows
%                     better abilities in structuture estimation than first
%                     estimate full (with families) heterogeneous binary
%                     HAC and then collapse it, due to a lot of biasing by
%                     the restrictions given by the SNC.
% Reestimator       - Applies if PreCollapse == true. The possible values 
%                     are either 'Ktauavg' or 'taumin', each corresponding to
%                     an approach to re-estimation of HAC's parameters
%                     after collapsing, see [Gorecki et al., 2017] for
%                     details.
% CollapsedArray    - Allows for supplying a collapsed array of '?'-homogenous
%                     HACs to speed up the estimating process.
% MinDistanceArray  - With CollapsedArray, provide also MinDistanceArray
%                     obtained by the method collapse.
% nForks            - An iteger from {1, ..., d-1} or 'unknown'. While
%                     collapsing, nForks determines, which collapsed HAC
%                     will be chosen as the final estimate. If an integer
%                     is provided, the collapsed HAC will have the number
%                     of forks equal to this integer. If set to 'unknown',
%                     the number of forks is determined using the
%                     function findjump, see Section 6.1 in [Gorecki et al., 2017].
%                     
%
% Example with the default settings of the optional parameters (omitting
% the ones that serve for supplying some precomputed quantities):
% [estimate, fitLog] = HACopulafit(U, families, ...
%                       'HACEstimator', 'pairwise',... 
%                       'ThetaEstimator', 'invtau', ...
%                       'ThetaEstimator2', 'invtau', ...
%                       'g_1', 'average', 'g_2', @(t)mean(t), 'GOF', 'R', ...
%                       'PreCollapse', true, 'Reestimator', 'Ktauavg', ...
%                       'nForks', 'unknown', 'Attitude', 'optimistic', ...
%                       'CheckData', 'on');
%
% Output (obligatory):
% HACObject     - An instance of the HACopula class.
% Output (optional):
% fitLog        - A string that describes the estimation process.
%
% NOTE: if a NaN or Inf is generated in the diagonal transformation, it is
%         replaced using the nanapprox function.
% 
%
% References:
% [Gorecki et al., 2016a] Gorecki, J., Hofert, M., and Holena, M. (2016). An 
%     approach to structure determination and estimation of hierarchical
%     Archimedean copulas and its application to bayesian classication.
%     Journal of Intelligent Information Systems, pages 21-59.
%
% [Gorecki et al., 2017] On Structure, Family and Parameter Estimation
%     of Hierarchical Archimedean copulas. Journal of Statistical Computation 
%     and Simulation, 87(17), 3261-3324
%
%
% Copyright 2018 Jan Gorecki

% -------------------------------------------------------------------------


%% the default setting of the parameters
hacEstimatorRatio = 0; % 'pairwise'
thetaEstimatorPairwise = 'invtau';
thetaEstimatorDiagonal = 'invtau';
g1 = @(t)mean(t);
g1name = 'average';
g2 = @(t)mean(t);
gof = 'R';
preCollapse = true; % do pre-collapsing?
reEstimType = 'Ktauavg';
nForks = 'unknown';
attitude = 'optimistic'; 
checkData = true;
% the following default parameters are the ones that can be somehow
% pre-computed (to increase efficinecy of the estimation process)
preCollapsedHAC = []; 
emp2copulas = {};
collapsedArray = {};
minDistanceArray = [];


%% inputs check

narginchk(2,34);

% switch on logging if the second output argument is specified
if nargout == 2
    logging = true;
elseif nargout <= 1
    logging = false;
end


% check for additional parameters
if mod(size(varargin,2),2) == 1
    error('HACopula:HACopulafit:BadInputs', 'HACopulafit: There must be an even number of the additional paramters');
end

if size(varargin,2) >= 2
    % there are some additional parameters
    parNames = lower(varargin(1:2:size(varargin,2)-1));
    parValues = varargin(2:2:size(varargin,2));
    
    % check of the parameter names are allowed
    for i = 1:size(parNames,2)
        try
            validatestring(parNames{i},{'ThetaEstimator', 'ThetaEstimator2', ...
                'GOF', 'Emp2copulas', 'Logging', 'g_1', 'g_2', 'CheckData', 'HACEstimator', ...
                'KendallMatrix', 'Attitude', 'Reestimator', 'nForks', 'CollapsedArray', 'MinDistanceArray', 'PreCollapse', 'PreCollapsedHAC'});
        catch
            error('HACopula:HACopulafit:BadInputs', ['The input, ''' parNames{i} ''', did not match any of the valid parameter names '...
                '(ThetaEstimator, ThetaEstimator2, GOF, Emp2copulas, g_1, g_2, CheckData'...
                ', HACEstimator, KendallMatrix, Attitude, Reestimator, nForks, CollapsedArray, MinDistanceArray, PreCollapse, PreCollapsedHAC).']);
        end
    end
    
    % check for ThetaEstimator parameter
    iParameter = find(strcmp(parNames, lower('ThetaEstimator')));
    % is ThetaEstimator a parameter ?
    if size(iParameter,2) > 0 %
        if size(iParameter,2) > 1
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: ThetaEstimator is a repeating parameter.')
        end
        thetaEstimatorPairwise = parValues{iParameter};
        if ~(strcmp(thetaEstimatorPairwise, 'invtau') || strcmp(thetaEstimatorPairwise, 'invtau2') || ...
                strcmp(thetaEstimatorPairwise, 'mle') || strcmp(thetaEstimatorPairwise, 'mle2'))
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: the value corresponding to ThetaEstimator must be ''invtau'', ''invtau2'', ''mle'' or ''mle2''.');
        end
        % set the diagonal theta estimator to the same value (which can be changed by setting ThetaEstimator2)
        thetaEstimatorDiagonal = parValues{iParameter};
    end
    
    % check for ThetaEstimator2 parameter (necessary in the case that combination of estimators are used)
    iParameter = find(strcmp(parNames, lower('ThetaEstimator2')));
    % is ThetaEstimator2 a parameter ?
    if size(iParameter,2) > 0 %
        if size(iParameter,2) > 1
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: ThetaEstimator2 is a repeating parameter.')
        end
        thetaEstimatorDiagonal = parValues{iParameter};
        if ~(strcmp(thetaEstimatorDiagonal, 'invtau') || strcmp(thetaEstimatorDiagonal, 'invtau2') || ...
                strcmp(thetaEstimatorDiagonal, 'mle') || strcmp(thetaEstimatorPairwise, 'mle2'))
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: the value corresponding to ThetaEstimator2 must be ''invtau'', ''invtau2'', ''mle'' or ''mle2''.');
        end
    end
    
    
    % check for GOF parameter
    iParameter = find(strcmp(parNames, lower('GOF')));
    % is GOF a parameter ?
    if size(iParameter,2) > 0 %
        if size(iParameter,2) > 1
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: GOF is a repeating parameter.')
        end
        gof = parValues{iParameter};
        if ~(strcmp(gof, 'E') || strcmp(gof, 'K') || ...
                strcmp(gof, 'R'))
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: the value corresponding to GOF must be ''E'', ''K'' or ''R''.');
        end
    end
    
    % check for Emp2copulas parameter
    iParameter = find(strcmp(parNames, lower('Emp2copulas')));
    % is Emp2copulas a parameter ?
    if size(iParameter,2) > 0 %
        if size(iParameter,2) > 1
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: Emp2copulas is a repeating parameter.')
        end
        emp2copulas = parValues{iParameter};
        if ~iscell(emp2copulas)
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: the value corresponding to Emp2copulas must a cell returned by computeallemp2copulas.');
        end
    end
    
    % check for g_1 parameter
    iParameter = find(strcmp(parNames, lower('g_1')));
    % is g_1 a parameter ?
    if size(iParameter,2) > 0 %
        if size(iParameter,2) > 1
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: g_1 is a repeating parameter.')
        end
        g1name = parValues{iParameter};
        % check if the parameter is max, min or average
        if strcmp(g1name,'min')
            g1 = @(t)min(t);
        elseif strcmp(g1name,'max')
            g1 = @(t)max(t);
        elseif strcmp(g1name,'average');
            g1 = @(t)mean(t);
        else
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: The value of g_1 is not a supported aggregation function. Use ''max'', ''min'' or ''average''');
        end
    end
    
    % check for g_2 parameter
    iParameter = find(strcmp(parNames, lower('g_2')));
    % is g_2 a parameter ?
    if size(iParameter,2) > 0 %
        if size(iParameter,2) > 1
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: g_2 is a repeating parameter.')
        end
        g2 = parValues{iParameter};
        % check if g_2 is a [0, Inf)-aggregation function for some
        % values
        if ~isnumeric(g2(0.5))
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: g_2(0.5) is not a number.')
        end
        if g2([0.25 0.75]) ~= g2([0.75 0.25])
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: g_2(0.25, 0.75) ~= g_2(0.75, 0.25), i.e., g_2 is not invariant to a permutation of arguments.')
        end
        if g2([0.5 0.5]) ~= g2(0.5)
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: g_2(0.5, 0.5) ~= g_2(0.5), i.e., g_2 is not a supported aggregation function.')
        end
    end
    
    
    % check for CheckData parameter
    iParameter = find(strcmp(parNames, lower('CheckData')));
    % is CheckData a parameter ?
    if size(iParameter,2) > 0 %
        if size(iParameter,2) > 1
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: CheckData is a repeating parameter.')
        end
        checkData = parValues{iParameter};
        if ~islogical(checkData) % checkData must be a logical value or
            % on or off
            if ~(strcmp(checkData, 'on') || strcmp(checkData, 'off'))
                error('HACopula:HACopulafit:BadInputs', 'HACopulafit: the value corresponding to CheckData must be ''on'', ''off'' or a logical value.');
            end
            checkData = strcmp(checkData, 'on');
        end
    end
    
    
    % check for HACEstimator parameter 
    iParameter = find(strcmp(parNames, lower('HACEstimator')));
    % is HACEstimator a parameter ?
    if size(iParameter,2) > 0 %
        if size(iParameter,2) > 1
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: HACEstimator is a repeating parameter.')
        end
        hacEstimatorRatio = parValues{iParameter}; % ratio between pairwise and diagonal HAC estimator
        if ~isnumeric(hacEstimatorRatio)
            % could be 'pairwise' or 'diagonal'
            if ~(strcmp(hacEstimatorRatio, 'pairwise') || strcmp(hacEstimatorRatio, 'diagonal'))
                error('HACopula:HACopulafit:BadInputs', 'HACopulafit: the string value corresponding to HACEstimator must be ''pairwise'' or ''diagonal''.');
            end
            if strcmp(hacEstimatorRatio, 'pairwise')
                hacEstimatorRatio = 0;
            else % hacEstimator = 'diagonal'
                hacEstimatorRatio = 1;
            end
        else % HACEstimator is numeric
            if (hacEstimatorRatio < 0) || (hacEstimatorRatio > 1)
                error('HACopula:HACopulafit:BadInputs', 'HACopulafit: the numeric value corresponding to HACEstimator must from the interval [0, 1].');
            end
        end
    end
    
    % check for KendallMatrix parameter
    iParameter = find(strcmp(parNames, lower('KendallMatrix')));
    % is KendallMatrix a parameter ?
    if size(iParameter,2) > 0 %
        if size(iParameter,2) > 1
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: KendallMatrix is a repeating parameter.')
        end
        precomputedKendallMatrix = parValues{iParameter};
    end
    
    % check for Attitude parameter
    iParameter = find(strcmp(parNames, lower('Attitude')));
    % is Attitude a parameter ?
    if size(iParameter,2) > 0 %
        if size(iParameter,2) > 1
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: Attitude is a repeating parameter.')
        end
        attitude = parValues{iParameter};
        if ~(strcmp(attitude, 'pessimistic') || strcmp(attitude, 'optimistic'))
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: the value corresponding to Attitude must be ''pessimistic'' or ''optimistic''.');
        end
    end
    
    % check for Reestimator parameter
    iParameter = find(strcmp(parNames, lower('Reestimator')));
    % is Reestimator a parameter ?
    if size(iParameter,2) > 0 %
        if size(iParameter,2) > 1
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: Reestimator is a repeating parameter.')
        end
        reEstimType = parValues{iParameter};
        if ~(strcmp(reEstimType, 'Ktauavg') || strcmp(reEstimType, 'taumin'))
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: the value corresponding to Attitude must be ''Ktauavg'' or ''taumin''.');
        end
    end
    
    % check for nForks parameter
    iParameter = find(strcmp(parNames, lower('nForks')));
    % is Reestimator a parameter ?
    if size(iParameter,2) > 0 %
        if size(iParameter,2) > 1
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: nForks is a repeating parameter.')
        end
        nForks = parValues{iParameter};
        if (isnumeric(nForks) && floor(nForks) ~= nForks) || (~isnumeric(nForks) && ~strcmp(nForks, 'unknown'))
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: the value corresponding to nForks must be an integer or ''unknown''.');
        end
    end
    
    % check for CollapsedArray parameter
    iParameter = find(strcmp(parNames, lower('CollapsedArray')));
    % is CollapsedArray a parameter ?
    if size(iParameter,2) > 0 %
        if size(iParameter,2) > 1
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: CollapsedArray is a repeating parameter.')
        end
        collapsedArray = parValues{iParameter};
        if ~isa(collapsedArray{1},'HACopula')
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: the value corresponding to CollapsedArray must be a cell array of HACopula objects.');
        end
    end
    
    % check for MinDistanceArray parameter
    iParameter = find(strcmp(parNames, lower('MinDistanceArray')));
    % is MinDistanceArray a parameter ?
    if size(iParameter,2) > 0 %
        if size(iParameter,2) > 1
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: MinDistanceArray is a repeating parameter.')
        end
        minDistanceArray = parValues{iParameter};
        if ~isnumeric(minDistanceArray)
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: the value corresponding to MinDistanceArray must be a numeric array.');
        else
            if ~isnumeric(nForks) && strcmp(nForks, 'unknown') && size(minDistanceArray,2)~=length(collapsedArray)
                error('HACopula:HACopulafit:BadInputs', 'HACopulafit: MinDistanceArray must be of the same size as the CollapsedArray.');
            end
        end
    end
    
    % check for PreCollapse parameter
    iParameter = find(strcmp(parNames, lower('PreCollapse')));
    % is PreCollapse a parameter ?
    if size(iParameter,2) > 0 %
        if size(iParameter,2) > 1
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: PreCollapse is a repeating parameter.')
        end
        preCollapse = parValues{iParameter};
        if ~islogical(preCollapse)
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: the value corresponding to PreCollapse must be a logical value (true or false).');
        end
    end
    
    % check for PreCollapsedHAC parameter
    iParameter = find(strcmp(parNames, lower('PreCollapsedHAC')));
    % is PreCollapsedHAC a parameter ?
    if size(iParameter,2) > 0 %
        if size(iParameter,2) > 1
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: PreCollapsedHAC is a repeating parameter.')
        end
        preCollapsedHAC = parValues{iParameter};
        if ~isa(preCollapsedHAC, 'HACopula')
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: the value corresponding to PreCollapsedHAC must be an instance of the HACopula class.');
        end
    end
    
    
end % more than two input parameters

% families check
ALL_FAMILIES = {'A', 'C', 'F', 'G', 'J', '12', '14', '19', '20', '?'}; % supported families
for i = 1:length(families)
    if sum(strcmp(families{i}, ALL_FAMILIES)) == 0
        error('HACopula:HACopulafit:BadInputs', ['HACopulafit: The provided family ''' families{i} ''' is not supported. Choose families from {' strjoin(ALL_FAMILIES, ', ') '}.']);
    end
end
if sum(strcmp('?', families)) >= 1 && length(families) > 1
    error('HACopula:HACopulafit:BadInputs', 'HACopulafit: The arbitrary family ''?'' cannot be combined with other families.');
end

% pre-collapsing check
if preCollapse && hacEstimatorRatio > 0
    error('HACopula:HACopulafit:BadInputs', ['HACopulafit: Using the diagonal estimation (i.e., ''HACEstimator'' == ''diagonal'' or ''HACEstimator'' > 0),' ... 
        ' pre-collapsing cannot be performed as pre-collapsing requires no assumptions on the underlying families (on contrary to the diagonal estimation).' ...
        ' Either set the parameter ''PreCollapse'' to false, or ''HACEstimator'' to 0 (''pairwise'').']);
end


% perform basic data checks
if checkData && ~iscopuladata(U)
    warning('HACopula:notCopulaData', 'HACopulafit: at least one column of U was rejected to be standard uniform according to KS test at the 5%% significance level, i.e., U might not be a sample from a copula.');
end

    
% for the heterogeneus case, compute empirical bivariate copulas if not provided
if (size(families, 2) > 1) && isempty(emp2copulas)
    emp2copulas = computeallemp2copulas(U);
end

%% -------------------------------------------------------------------------
% Initialization steps:

[n, d] = size(U);

%initialize nesting
N = getN0(families);
N1 = cell(1,d);
for i = 1:d
    N1{i} = N;
end

%initialize the descent leaves (\downarrow(i), i = 1, ..., d)
descLeaves = cell(1,d);
for i = 1:d
    descLeaves{i} = i;
end

% initialize the set I
I = 1:d;

% set Kendall matrix
if exist('precomputedKendallMatrix', 'var')  % NOTE: the term eps(1) appears two times to match compatibility with Octave
    if ~(ismatrix(precomputedKendallMatrix) && size(precomputedKendallMatrix, 1) == d && ...
        size(precomputedKendallMatrix, 2) == d && sum(diag(precomputedKendallMatrix) - 1) <= d*eps(1) && ...
        sum(sum( (precomputedKendallMatrix >= -1) .* (precomputedKendallMatrix <= (1+eps(1))  ) )) == d^2 )
            error('HACopula:HACopulafit:BadInputs', 'HACopulafit: the value corresponding to KendallMatrix is not a Kendall correlation matrix.');
    end
    K = precomputedKendallMatrix;
else
    %compute the sample version of the Kendall correlation matrix
    K = kendallTauMatrix(U);
end

% check for negative pairwise correlation
% NOTE: if the average of the negative pairwise correlation coefficients in
% K is lower than MIN_AVERAGE_TAU, a warning is shown
if checkData
   MIN_AVERAGE_TAU = -0.1; 
   negK = K < 0; % identify negative coefficients in K
   if sum(sum(negK)) > 0
       avgNegK = sum(sum(K(negK)))/sum(sum(negK));  % compute the average of the negative coefficients
       if avgNegK < MIN_AVERAGE_TAU
           warning('HACopula:HACopulafit', ['HACopulafit: The average of the negative coefficient in the Kendall''s matrix is lower that ' ...
               num2str(MIN_AVERAGE_TAU) ' . Recall that HACs are not able to model negative pairwise correlation.' ...
               ' Try to flip (X := - X, using the function findvars2flip) some of the variables in order to reduce the negative pairwise correlation and possibly better fit of the estimate.']);
       end
   end
end
            
    
%--------------------------------------------------------------------------
% STARTING ESTIMATION PROCESS
%--------------------------------------------------------------------------
%start logging
if logging
    fitLog = ['***** HACopulafit: Start... *****' char(10)];
end
isTauMatInLog = false;

%start the HAC estimation process
if ~isempty(preCollapsedHAC)
    colEstimate = preCollapsedHAC;
    forksToEst = length(colEstimate.Forks);
else
    % pre-collapsing of the structure
    doCollapsing = (preCollapse && hacEstimatorRatio == 0); % pairwise hacEstimator && the heterogeneous case
    if doCollapsing
        if isempty(collapsedArray)
            % estimate a binary ?-homogenous HAC
            HOMO_FAM = {'?'}; % an arbitrary family that serves just for precollapsing (allows for all taus in [-1, 1])
            if logging
                [binEstimate, fitLogHomo] = HACopulafit(U, HOMO_FAM, 'HACEstimator', 0,...
                    'g_1', g1name, 'Attitude', attitude, 'CheckData', false, 'KendallMatrix', K,...
                    'PreCollapse', false);
            else
                binEstimate = HACopulafit(U, HOMO_FAM, 'HACEstimator', 0,...
                    'g_1', g1name, 'Attitude', attitude, 'CheckData', false, 'KendallMatrix', K,...
                    'PreCollapse', false);
            end
            if logging
                fitLog = [fitLog '***** Pre-collapsing: Start... *****' char(10) ...
                          'Estimating the binary structure (no assumptions on the families)...' char(10) ...
                    fitLogHomo char(10)];
                isTauMatInLog = true;
            end            
            [colEstArray, minDistArray] = collapse(binEstimate, 'invtau', K, U, @(t)mean(t), 'optimistic', reEstimType, false);
            if logging
                fitLog = [fitLog char(10) 'Generating a sequence of ' num2str(length(colEstArray)) ' (collapsed) structures from the binary one.' char(10)];
            end            
        else
            colEstArray = collapsedArray;
            minDistArray = minDistanceArray;
        end
        % get(estimate) the number of forks and get the corresponding collapsed
        % homogeneous HAC
        if isnumeric(nForks) % k-known assumption
            if nForks <= d-1 && nForks >= 1
                colEstimate = colEstArray{min(d-nForks, length(colEstArray))}; % the min is used due to the possible case that considerFamilies = true in the collapse function
            else
                error('HACopula:HACopulafit:BadInputs', 'HACopulafit: the parameter nForks must be from {1, ..., d-1}');
            end
            if logging
                fitLog = [fitLog 'Taking the collapsed structure with ' num2str(nForks) ' forks, ' char(10) ...
                          'following the input ''nForks'', i.e., ' char(10)];
            end
        else % k-unknown assumption
            iJump = findjump(minDistArray); % using a heuristic to estimate the true number of the forks
            colEstimate = colEstArray{min(iJump, length(colEstArray))}; % if not collapsed to iJump forks, then the one with the lowest number of forks
            if logging
                fitLog = [fitLog 'Taking the collapsed structure with ' num2str(d-min(iJump, length(colEstArray))) ' forks, ' char(10) ...
                          'following the estimated number of forks given by findjump, i.e.,' char(10)];
            end
        end
        forksToEst = length(colEstimate.Forks); % estimate collapsed forks only
        if logging 
            fitLog = [fitLog 'instead of the binary structure given by' char(10)];
            for iForks = 1:d-1
                frk = getforkwithtauindex(colEstArray{1}, d+iForks);
                fitLog = [fitLog '\downarrow(' num2str(d+iForks) ') = {' num2str(sort(frk.Leaves)) '}' char(10)];
            end
            fitLog = [fitLog 'taking the non-binary one given by' char(10)];
            for iForks = 1:forksToEst
                frk = getforkwithtauindex(colEstimate, d+iForks);
                fitLog = [fitLog '\downarrow(' num2str(d+iForks) ') = {' num2str(sort(frk.Leaves)) '}' char(10)];
            end
            fitLog = [fitLog 'Note that \downarrow(i) = {i} for i = 1, ..., ' num2str(d) '.' char(10)];
            fitLog = [fitLog '***** Pre-collapsing: Done. *****' char(10) char(10) 'Estimating the families and parameters...'];
        end
    else
        forksToEst = d - 1; % estimate d - 1 forks, i.e., a binary HAC
    end
end


if logging && ~isTauMatInLog
    fitLog = [fitLog '(tau^n_{ij}):' char(10)];
    for i = 1:size(K,1)
        for j = 1:size(K,2)
            if i < j
                fitLog = [fitLog sprintf('%2.4f',K(i,j)) ' '];
            else
                fitLog = [fitLog '   .   '];
            end
        end
        fitLog = [fitLog char(10)];
    end
    fitLog = fitLog(1:end-1);
end


HACCellStruc = cell(1, 2*d-1);
for i = 1:d
    HACCellStruc{i} = i;
end

for k = 1:forksToEst
    if logging
        fitLog = [fitLog char(10) '--------------------------------------------------------------------' char(10)];
        fitLog = [fitLog 'k = ' num2str(k) ' *** Estimating \lambda(' num2str(k+d) '):' char(10)];
    end
    
    if ~isempty(preCollapsedHAC) || doCollapsing % we already have a pre-collapsed HAC
        fork = getforkwithtauindex(colEstimate, d+k);
        i = zeros(1, length(fork.Child));
        for iChild = 1:length(fork.Child)
            if isa(fork.Child{iChild}, 'HACopula')
                i(iChild) = fork.Child{iChild}.TauOrdering;
            else
                i(iChild) = fork.Child{iChild};
            end
        end
        maxSimilarity = fork.Tau;
    else % estimate a binary structured HAC
        %find maximal similarity in order to estimate the HAC structure
        maxSimilarity = -Inf;
        for ii = I              % ii is \tilde{i} in the paper
            for jj = I          % jj is \tilde{j} in the paper
                if (ii < jj)
                    if hacEstimatorRatio == 0 % pairwise hacEstimator
                        % aggregate the values of the Kendall correlation matrix
                        subK = K(descLeaves{ii}, descLeaves{jj});
                        vecSubK = reshape(subK, [1, length(descLeaves{ii})*length(descLeaves{jj})]);
                        similarity = g1(vecSubK);
                    elseif hacEstimatorRatio == 1 % the diagonal hacEstimator
                        similarity = K(ii, jj); % the K matrix is enlarged (in later steps) by one in each step
                    else % hacEstimator in (0, 1)
                        subK = K(descLeaves{ii}, descLeaves{jj});
                        vecSubK = reshape(subK, [1, length(descLeaves{ii})*length(descLeaves{jj})]);
                        similarityPairwise = g1(vecSubK);
                        similarityDiagonal = K(ii, jj); % the K matrix is enlarged (in later steps) by one in each step
                        similarity = (1-hacEstimatorRatio) * similarityPairwise + hacEstimatorRatio * similarityDiagonal; % weighted average
                    end
                    if (similarity > maxSimilarity)
                        maxSimilarity = similarity;
                        ii_max = ii; jj_max = jj;
                    end
                end
            end
        end
        %(i, j) from the paper
        i = [ii_max jj_max];
    end
    
    %nesting
    N = N1{i(1)};
    for iii = 2:size(i, 2)
        N = intersectNs(N, N1{i(iii)});
    end
    
    
    if logging
        fitLog = [fitLog 'I = [' num2str(I) ']' char(10)];
        %fitLog = [fitLog 'argmax g1(K_ii,jj) (i, j) = (' num2str(i) ', ' num2str(j) ')' char(10)];
        %fitLog = [fitLog '\downarrow(' num2str(i) ') = [' num2str(descLeaves{i}) '] \downarrow(' num2str(j) ') = [' num2str(descLeaves{j}) ']' char(10)];
        fitLog = [fitLog 'g_1(K_{'];
        sorI = sort(i);
        for dwnI = 1:size(i,2)
            fitLog = [fitLog '\downarrow(' num2str(sorI(dwnI)) ')'];
            if dwnI < size(i,2)
                fitLog = [fitLog ', '];
            end
        end
        fitLog = [fitLog '}) = ' num2str(maxSimilarity) char(10)];
        
        try
        if ~strcmp(N{1}{1}, '?') % omit for pre-collapsing
            %fitLog = [fitLog char(10)];
            fitLog = [fitLog 'Admissible families + ranges: ' char(10) '{' ];
            for iN = 1:length(N)
                if isinf(N{iN}{2}(2))
                    prnth = ')';
                else
                    prnth = ']';
                end
                if N{iN}{2}(1) == eps(0)
                    lBound = 'eps(0)';
                else
                    lBound = num2str(N{iN}{2}(1));
                end
                fitLog = [fitLog '(' N{iN}{1} ', [' lBound ', ' ...
                         num2str(N{iN}{2}(2)) prnth ')'];
                if iN < length(N)
                    fitLog = [fitLog ', '];
                end
            end
            fitLog = [fitLog '}' char(10)];
            fitLog = [fitLog 'Theta estimation:' char(10)];
        end
        catch
            bb =0;
        end
    end
    
    % estimate the parameters for all families in N
    Sn = zeros(1, length(N));
    theta = zeros(1, length(N)); % estimate as many thetas as families in N
    nDisqualified = 0; % counts how many families are disqualified from the estimation process
    for l = 1:length(N)
        a = N{l}{1};    % family and
        r = N{l}{2};    % its parameter range satisfying the s.n.c.
        
        % estimate the 2-AC parameter
        thetaAct = fitparameter(descLeaves(i), a, ...
            hacEstimatorRatio, thetaEstimatorPairwise, thetaEstimatorDiagonal, ...
            K, U, g1, i(1), i(2), maxSimilarity, r, attitude);
        
                
        if isnan(thetaAct) && logging
            fitLog = [fitLog 'no estimates available; family ' a ' is disqualified from the process' char(10)];
        end
        
        if logging && ~strcmp(N{1}{1}, '?') % omit for pre-collapsing
            fitLog = [fitLog 'family = ' a '   theta = ' num2str(thetaAct) char(10)];
        end
        
        % use GOF to find the best fitting family
        if (~isnan(thetaAct)) % if the family is not disqualified from the process
            thetaTrimmed = min(max(r(1),thetaAct),r(2));
            isThetaTrimmed = false;
            if (thetaAct ~= thetaTrimmed)
                isThetaTrimmed = true;
                if logging
                    if strcmp(attitude, 'optimistic')
                        fitLog = [fitLog 'trimming theta to: ' num2str(thetaTrimmed) ' , i.e.,' char(10)];
                        fitLog = [fitLog 'family = ' a '   theta = ' num2str(thetaTrimmed) char(10)];
                    else %attitude == pessimistic
                        fitLog = [fitLog '(' a ', ' num2str(thetaAct) ') is not admissible for \lambda(' num2str(d+k) '). Disqualifying...' char(10)];
                    end
                end
                thetaAct = thetaTrimmed;
            end
            
            if isThetaTrimmed && ~strcmp(attitude, 'optimistic') % i.e., pessimistic
                % diqualified from the process
                nDisqualified = nDisqualified + 1;
                Sn(l) = Inf;
            else
                theta(l) = thetaAct;
                
                %evaluate GOF of \psi^(a, theta)
                if (length(N) >= 2)
                    g2SnPair = 0;
                    g2SnDiag = 0;
                    if hacEstimatorRatio < 1 % pairwise or combined
                        aggSnMatPair = [];
                        for iii = 1:size(i, 2)-1
                            for jjj = iii+1:size(i, 2)
                                SnMat = zeros(max(descLeaves{i(iii)}), max(descLeaves{i(jjj)}));
                                for ii = descLeaves{i(iii)}
                                    for jj = descLeaves{i(jjj)}
                                        SnMat(ii,jj) = gof2Sn(U(:, [ii jj]), a, theta(l), gof, emp2copulas{min(ii,jj),max(ii, jj)});
                                    end
                                end
                                aggSnMatPair = [aggSnMatPair reshape(SnMat(descLeaves{i(iii)}, descLeaves{i(jjj)}), [1, length(descLeaves{i(iii)})*length(descLeaves{i(jjj)})])];
                            end
                        end
                        g2SnPair = g2(aggSnMatPair);
                    end
                    if hacEstimatorRatio > 0 % diagonal or combined
                        aggSnMatDiag = [];
                        for iii = 1:size(i, 2)-1
                            for jjj = iii+1:size(i, 2)
                                indI = min(i(iii), i(jjj));
                                indJ = max(i(iii), i(jjj));
                                if indJ > size(emp2copulas,2) || isempty(emp2copulas{indI, indJ}) % add new emp2copulas if needed
                                    emp2C = computeallemp2copulas(U(:, [indI indJ]));
                                    emp2copulas(indI,indJ) = emp2C(1,2);
                                    emp2copulas{indJ,indI} = [];
                                end
                                SnScalar = gof2Sn(U(:, [indI indJ]), a, theta(l), gof, emp2copulas{indI,indJ});
                                aggSnMatDiag = [aggSnMatDiag SnScalar];
                            end
                        end
                        g2SnDiag = g2(aggSnMatDiag);
                    end
                    
                    %use the aggregation function g2 to aggregate the
                    %evaluations for the bivariate margins
                    Sn(l) = (1-hacEstimatorRatio) * g2SnPair  + hacEstimatorRatio * g2SnDiag; % weighted average
                    % NOTE: for binary HACs, diagonal is independent on g2
                else
                    %if #N = 1, no GOF is necessary
                    Sn(l) = 0;
                end
            end
        else % diqualified from the process
            nDisqualified = nDisqualified + 1;
            Sn(l) = Inf;
        end
    end
    
    if nDisqualified == length(N)
        warning('HACopula:HACopulafit', 'HACopulafit: all families disqualified from the estimation process. Returning no HAC (i.e., the empty matrix []). Use another (wider set of) families for the data.');
        if logging
            fitLog = [fitLog char(10) '***** HACopulafit: Stop. *****' char(10) ...
                'All families have been disqualified from the estimation process.' char(10)...  '
                'Rejecting to return a HAC. Returning the empty matrix [].' char(10) ...
                'Try to use another (wider set of) families or the optimistic attitude.'];
        end
        HACObject = [];
        return
    end
    
   
    % find the best fitting family
    [~, iBest] = min(Sn);
    family = N{iBest}{1};
    parameter = theta(iBest);
    if logging
        Nnames = cell(1, length(N));
        for iN = 1:length(N)
            Nnames{iN} = N{iN}{1};
        end
        SnChar = sprintf('%2.4f, ',Sn);
        SnChar = SnChar(1:end-2);
        if ~strcmp(N{1}{1}, '?') % omit for pre-collapsing
            fitLog = [fitLog 'Family estimation:' char(10)];
            fitLog = [fitLog 'S_n^{g_2} for the families (' strjoin(Nnames, ', ') ') is ('  SnChar ')' char(10)];
            fitLog = [fitLog 'Best-fitting \psi^(family, theta) = \psi^(' family ', ' num2str(parameter) ')' char(10)];
        end
    end
    
    %nesting
    N1{d + k} = intersectNs(N, computeN2(family, parameter, families));
    
    %exclude clustered
    I = setdiff(union(I, d + k),i);
    descLeaves{d+k} = [descLeaves{i}]; % do a union of children's leaves
    if logging
        fitLog = [fitLog 'Join leaves: \downarrow(' num2str(d+k) ') = [' num2str(sort(descLeaves{d+k})) ']'];
    end
    
    %build structure
    HACCellStruc{d+k} = {{family, parameter}, HACCellStruc{i}};
    
    % extend U and K for the diagonal HAC estimator
    if hacEstimatorRatio > 0 % combined or diagonal HAC estimator
        % construct the diagonal transformation
        psi = getgenerator(family, parameter, 0);
        psiInv = getgenerator(family, parameter, 1);
        delta = @(t) psi(size(i,2).* psiInv(t));  % the transformation
        U(:,d+k) = delta(max(U(:,i),[], 2)); % extend U
        % NaN check
        [U(:,d+k), nNaNs] = nanapprox(U(:,d+k), U(:,i));
        if nNaNs > 0
            warning('HACopula:NaN_detected', ['HACopulafit: ' num2str(nNaNs) ' NaNs detected in the diagonal transformation and replaced by their approximations.']);
        end
        % extend K with taus for the new vector U_{d+k}
        for ii = 1:d-k-1  % length(I) without the last one
            tau = kendallTauMatrix(U(:,[I(ii) d+k]));
            K(I(ii), d+k) = tau(1, 2);
            K(d+k, I(ii)) = tau(1, 2);
        end
        K(d+k, d+k) = 1;
    end
end

HACObject = HACopula(HACCellStruc{d+forksToEst});

if logging
    fitLog = [fitLog char(10) '--------------------------------------------------------------------' char(10) '***** HACopulafit: Stop. *****'];
end

end




