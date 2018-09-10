% The high-dimensional example
%
%
% Copyright 2018 Jan Gorecki

HACModel = getfullmodel(11, 10, 'C', 0.1, 0.08);
plot(HACModel);

if isoctave
    disp('Loading data...'); 
    load('Demos/highdimex_data.mat'); % load U generated in MATLAB
    disp('Loaded.'); 
else % MATLAB    
    disp('Sampling...'); 
    rng('default'); 
    rng(1); % set the seed
    tic;
    U = pobs(rnd(HACModel, 2000)); % do sampling
    toc
end

% NOTE: In Octave, the (standard) computation of the Kendall matrix requires a lot of 
% memory and we have experienced several times a system (Windows) hangup
% during the computation. However, after closing all other applications, 
% the computation has finished successfully.
tic; disp('Slow (O(n^2)) computation of Kendall''s matrix...'); 
KSlow = kendallTauMatrix(U);
toc

% compile fastKendallTau.cpp to a MATLAB executable (.mex) file
if isoctave
     mkoctfile --mex fastKendallTau.cpp
else % MATLAB
    mex fastKendallTau.cpp
end

tic; disp('Fast (O(n*log(n))) computation of Kendall''s matrix...'); 
K = kendallTauMatrix(U);
toc

disp('Are the results the same?');
sum(sum(K == KSlow)) == numel(K)

tic; disp('Estimating (1)...');
fit1Avg = HACopulafit(U, {'C'}, 'Reestimator', 'Ktauavg', 'KendallMatrix', K);
toc

tic; disp('Estimating (2)...');
fit2Min = HACopulafit(U, {'C'}, 'Reestimator', 'taumin', 'KendallMatrix', K);
toc

% plot(fit1Avg);
% plot(fit2Min);


%% Evaluation

tic; disp('goodness-of-fit');
[gofdSnE(fit1Avg, U) ...
 gofdSnE(fit2Min, U)]
toc

tic; disp('kendall (HAC vs sample)');
[distance(fit1Avg, K) ...
 distance(fit2Min, K)]
toc

DISTANCE_TYPE = {'kendall', 'lower-tail'};
for i = 1:2
    tic; disp([DISTANCE_TYPE{i} ' (HAC vs HAC)']);
    [distance(fit1Avg, HACModel, DISTANCE_TYPE{i}) ...
     distance(fit2Min, HACModel, DISTANCE_TYPE{i})]
    toc
end

tic; disp('structures match ratio...')
[~, fit1AvgRatio] = comparestructures(HACModel, fit1Avg);
[~, fit2MinRatio] = comparestructures(HACModel, fit2Min);
[fit1AvgRatio fit2MinRatio]
toc

tic; disp('evaluating CDF at (0.5, ..., 0.5)...')
[cdf(fit1Avg, 0.5 * ones(1, HACModel.Dim)) ...
cdf(fit2Min, 0.5 * ones(1, HACModel.Dim))]
toc