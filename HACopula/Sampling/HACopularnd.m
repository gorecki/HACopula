function [U, I] = HACopularnd(HACModel, n, varargin)
% HACOPULARND - Sampling from a hierarchical Archimedean copula.
%
% Purpose:
% The function generates a sample of *n* observations from the HAC
% described by an instance of the HACopula class *HACModel*. All sampling
% algorithms are based on the approaches from [Hofert, 2010] and [Hofert,
% 2012].
% 
% Usage:
% U = HACopularnd(HACModel, n)
% Note that the second output argument and varargin are used just for inner
% purposes of the function.
%
% Inputs:
% HACModel  - An instance of the HACopula class.
% n         - The number of observations (rows) to be generated.
%
% Output:
% U         - An (n * d)-matrix of observations from the d-HAC given by HACModel.
%
% NOTE: 
% 1) (CAUTION) Sampling a HAC with parameters close to tau = 1 could be
% numerically unstable and NaNs could be generated - in this case, solution
% could be reduction of dimension (removing strongly correlated
% variables).
% 2) Generators with tau < TAU_INDEPENDENCE are substituted by the
% independence generator (to get numerical stability).
% 3) If the parameters of the generators are close to their parameter
% ranges, they are slightly shifted more far from the bounds in order to get
% numerically stable results (see the function adjustparameters01 for more
% details).
% 4) If some NaNs (Infs) occur during the sampling process (but the ratio
% of the rows with a NaN is lower than NAN_INF_LIMIT), HACopularnd tries to
% execute itself several times (but MAX_REPEATS-times at maximum) in order to
% generate proper observations from [0, 1]^d.
%
% References:
% [Hofert, 2010] Hofert, M. (2010). Sampling Nested Archimedean Copulas
%     with Applications to CDO Pricing. Suedwestdeutscher Verlag fuer
%     Hochschulschriften.
%
% [Hofert, 2012] Hofert, M. (2012). A stochastic representation and
%     sampling algorithm for nested Archimedean copulas. Journal of Statistical
%     Computation and Simulation, 82(9):1239-1255.
%
%
% Copyright 2017 Jan Górecki

% -------------------------------------------------------------------------

TAU_INDEPENDENCE = 1e-15; % generators with tau lower than TAU_INDEPENDENCE 
                          % are substituted by the independence generator

% is it variable?
if ~isa(HACModel, 'HACopula')
    % if so, sample standard uniform distribution
    U = rand(n,1);
    I = HACModel;
    return;
end

if HACModel.Level == 1
    %Marshal-Olkin
    if HACModel.Tau >= TAU_INDEPENDENCE
        % if necessary, adjust the parameter in order to get numerically stable results
        theta0adj = adjustparameter0(HACModel.Family, HACModel.Parameter);
        % sample V0 ~ F0
        v0 = rndV0(HACModel.Family, theta0adj, n);
    else
        v0 = [];
    end
    U = cell(1, size(HACModel.Child, 2));
    I = cell(size(U));
    for iChild = 1:size(HACModel.Child, 2)
        [U{iChild}, I{iChild}] = HACopularnd(HACModel.Child{iChild}, n, v0); 
    end

else
    %McNeil
    parentFamily = HACModel.Parent.Family; 
    parentParameter = HACModel.Parent.Parameter;
    parentTau = HACModel.Parent.Tau;
    v0 = varargin{1};
    if parentTau < TAU_INDEPENDENCE
        if HACModel.Tau >= TAU_INDEPENDENCE
            % if necessary, adjust the parameter in order to get numerically stable results
            theta0adj = adjustparameter0(HACModel.Family, HACModel.Parameter);
            % sample V0 ~ F0
            v01 = rndV0(HACModel.Family, theta0adj, n);
        else
            v01 = [];
        end
    else
        % if necessary, adjust the parameters in order to get numerically stable results
        [theta0adj, theta1adj] = adjustparameters01(parentFamily, HACModel.Family, parentParameter, HACModel.Parameter);
        % sample V01 ~ F01
        v01 = rndV01(parentFamily, HACModel.Family, theta0adj, theta1adj, v0);
    end
    
    %recursion
    U = cell(1, size(HACModel.Child, 2));
    I = cell(size(U));
    for iChild = 1:size(HACModel.Child, 2)
        [U{iChild}, I{iChild}] = HACopularnd(HACModel.Child{iChild}, n, v01);
    end
end

% join matrices
U = cell2mat(U);
I = cell2mat(I);

if HACModel.Level == 1
    if HACModel.Tau >= TAU_INDEPENDENCE 
        %Marshal-Olkin
        gen0 = getgenerator(HACModel.Family, theta0adj, 0);
        for i = 1:size(U,2)
            U(:,i) = gen0(-log(U(:,i)) ./ v0);
        end
    end
    %permute data columns according to the HAC structure
    U(:,I) = U;
else
    if parentTau < TAU_INDEPENDENCE
        if HACModel.Tau >= TAU_INDEPENDENCE
            %Marshal-Olkin
            gen0 = getgenerator(HACModel.Family, theta0adj, 0);
            for i = 1:size(U,2)
                U(:,i) = gen0(-log(U(:,i)) ./ v01);
            end
        %else
            %do nothing - assumed independence
        end
    else
        %McNeil
        gen01 = getnestedgenerators(parentFamily, HACModel.Family, theta0adj, theta1adj);
        for i = 1:size(U,2)
            U(:,i) = gen01(-log(U(:,i)) ./ v01, v0);
        end
    end
end

% checks
NAN_INF_LIMIT = 0.99; % 1 = switch off this check
MAX_REPEATS = 50;
% check for NaNs and Infs and there are less then (or equal to) NAN_INF_LIMIT*n then 
% sample another observations and replace the NaNs and Infs. Otherwise
% throw an error.
% NOTE: even using the following correction of NaNs and the functions
% adjustParameters01 and adjustParameter0, it might happen than for some
% (extreme) combination of parameters in a HAC, this function (HACopularnd)
% will not return a valid sample. In this case, contact the author.

cont = true;
nRepeats = 0;
while cont
    infRows = any(isinf(U) | isnan(U),2);
    nInfRows = sum(infRows);
    if nInfRows > 0
        if nInfRows/n <= NAN_INF_LIMIT
            % sample additional rows (observations)
            if HACModel.Level == 1
                Xadd = HACopularnd(HACModel, nInfRows);
            else
                varargin{1} = max(v0(infRows), eps); % NOTE: if v0 = 0, set v0 = eps - prevents from NaNs
                Xadd = HACopularnd(HACModel, nInfRows, varargin{:});
            end
            % replace the rows with inf by the new ones
            U(infRows,:) = Xadd;
        else
            error(['HACopularnd: maximal limit of rows containing NaNs or Infs exceeded, i.e., > ' num2str(NAN_INF_LIMIT*100) '% .']);
        end
        cont = sum(any(isinf(U) | isnan(U),2)) >= 1;
        nRepeats = nRepeats + 1;
        if nRepeats > MAX_REPEATS
            error('HACopularnd: maximal number of repeating of replacing inf rows exceeded.');
        end
    else
        cont = false;
    end
end

end

function [theta0, theta1] = adjustparameters01(family0, family1, theta0, theta1)

% first adjust separetely the theta0 parameter
theta0 =  adjustparameter0(family0, theta0);
% and the theta1 parameter
theta1 =  adjustparameter0(family1, theta1);

% then adjust combination of the parameters

% adjust homogeneous sampling
% avoid sampling for theta0 == theta1 for a homogenneous pair in order to
% get numerically stable results
% NOTE: if all pairs of child-parent
% nodes from the same family with theta0 ~=~ theta1 are collapsed to a
% single node, this adjustment is not necessary
if (strcmp(family0, family1)) &&  ((theta1 - theta0) < 1e-15)
    theta1 = theta0 + 1e-15;
end

familyNames = [family0 family1];
% adjust heterogeneous sampling
% the constants are set experimentally according to sampling a heterogenous
% 3-HAC and 4-HAC
switch familyNames
    %case 'AC'
        % no adjustments needed
    case {'A19', 'C19'}
        theta1 = max(theta1, 1e-12);
    case 'A20'
        theta1 = max(theta1, 1 + 1e-1); 
        % NOTE: rndV01 may produce non-uniform (rejected if alpha = 5%)
        % data even if eps = 1e-1, e.g., for HACopula({{'A', 0.5}, 1, {{'A', 0.9}, 2, {{'20', 1}, 3, 4}}});
        % however: this data are very close to uniform data, hence
        % repeating sampling several times produces uniform data in few
        % repeats
    case {'C12','C14'}
        theta1 = max(theta1, 1 + eps);        
    case 'C20'
        if (theta1 - theta0) < 1e-15 
            theta1 = theta0 + 1e-15;
        end
end
end


function theta0 = adjustparameter0(family0, theta0)

switch family0
    case 'A'
        theta0 = min(theta0, 1 - 1e-12);
    case {'12','14'}
        theta0 = max(theta0, 1 + eps);
    case '19'
        theta0 = max(theta0, 1e-12);        
end
end