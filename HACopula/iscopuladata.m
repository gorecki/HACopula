function [isCopulaData, hArray, pArray] = iscopuladata(U, varargin)
%ISCOPULADATA - Check for copula data
% 
% The function iscopuladata checks:
% 1) if there are no NaN or Inf in U (if not, an error is thrown)
% 2) if all elements in U are from [0, 1] (if not, an error is thrown)
% 3) if the observations U are close to standard uniform distribution (for
% each margin) - using two-sample Kolmogorov-Smirnov test against a sample
% with a perfect standard uniform distribution
%
% Usage:    [isCopulaData, hArray, pArray] = iscopuladata(U);
%           [isCopulaData, hArray, pArray] = iscopuladata(U, alpha)
%
% Input:
%           U               - (A (n x d) dimensional vector of values lying 
%                             in [0,1] (the observations).
%           alpha           - a significance level from [0, 1] for the two-sample
%                             Kolmogorov-Smirnov test (0.05 by default)
% Output:
%           isCopulaData    - isCopulaData = 1 if there is no rejection for all margins
%                             isCopulaData = 0 otherwise
%           hArray          - vector of results of the two-sample
%                             Kolmogorov-Smirnov test for each margin (0
%                             non-rejected, 1 - rejected)
%           pArray          - vector of p-values for the two-sample
%                             Kolmogorov-Smirnov test for each margin 
%
%
% Copyright 2018 Jan Gorecki

narginchk(1,2);

if size(varargin, 2) == 1
    alpha = varargin{1};
    if ~(isnumeric(alpha) && (alpha >= 0) && (alpha <= 1))
        error('HACopula:BadInputs', 'The paramter ''alpha'' must be a number in [0, 1].')
    end
else
    alpha = 0.05;
end

[n,d] = size(U);
if n < 3
    error('HACopula:iscopuladata', 'iscopuladata: the sample is too small (too few rows, i.e., n < 3)');
end

% check if there is no NaN in U
if (sum(sum(isnan(U))) > 0) || sum(sum(isinf(U)))
    error('HACopula:BadInputs', 'iscopuladata: There are NaN or Inf in the input matrix U');
end
    

% check if the data lies in the unit hypercupe
if sum(sum(U < 0)) || sum(sum(U > 1))
    error('HACopula:BadInputs', 'iscopuladata: The data has to be strictly in the intervall [0,1]')
elseif sum(sum(U == 0)) || sum(sum(U == 1))
    warning('HACopula:iscopuladata', 'iscopuladata: Some of the data is on the bound of the unit cube. Not all methods are stable or some functions give back non-finite results, if some input (Copula-)data is exactly 0 or 1.')
end

if isoctave
    warning('HACopula:iscopuladata', ['iscopuladata (for Octave): the ''kstest2'' function used in iscopuladata' ...
        ' to test marginal uniformity of the data belongs to the statistics package from Octave' ...
        'Forge but has not yet been implemented (9.8.2018). Omitting this test...']);
    isCopulaData = 1;
    hArray = zeros(1, d);
    pArray = zeros(1, d);
else % MATLAB
    
    % check if the sample is close a copula distribution
    %
    % generate a sample with a perfect standard uniform distribution
    perfUniform = (0:1/(n-1):1)'; %e.g., [0 0.25 0.5 0.75 1] for n = 5
    
    hArray = zeros(1, d);
    pArray = zeros(1, d);
    for i = 1:d
        % compare i-th margin with the sample with a perfect standard uniform distribution
        % use two-sample Kolmogorov-Smirnov test
        [hArray(i), pArray(i)] = kstest2(U(:,i), perfUniform, 'Alpha',alpha);
    end
    
    % isCopulaData = 1 if there is no rejection for all margins
    % isCopulaData = 0 otherwise
    isCopulaData = (min(sum(hArray), 1) == 0);
    
end
end