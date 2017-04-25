% The high-dimensional example
%
%
% Copyright 2017 Jan Górecki

HACModel = gethomomodel(11, 10, 'C', 0.1, 0.08);
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

tic; disp('Kendall''s matrix...'); 
if isoctave
    K = kendall(U);
else % MATLAB
    K = corr(U, 'type', 'kendall');
end
toc

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

tic; disp('comparing the structures...')
[comparestructures(HACModel, fit1Avg) ...
comparestructures(HACModel, fit2Min)]
toc

tic; disp('evaluating at (0.5, ..., 0.5)...')
[evaluate(fit1Avg, 0.5 * ones(1, getdimension(HACModel))) ...
evaluate(fit2Min, 0.5 * ones(1, getdimension(HACModel)))]
toc