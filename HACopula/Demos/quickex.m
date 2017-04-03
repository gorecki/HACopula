% The quick example 
%
%
% Copyright 2017 Jan Górecki


%% Construct & plot a HAC model
LAM8 = {'12', tau2theta('12', 0.8)}
LAM9 = {'19', tau2theta('19', 0.7)}
LAM10 = {'12', tau2theta('12', 0.5)}
LAM11 = {'C', tau2theta('C', 0.2)}

HACModel = HACopula({LAM11, {LAM9, 2, 5, 6}, {LAM10, 1, {LAM8, 3, 4, 7}}});

plot(HACModel);

%% Computing probabilities involving a HAC
disp('HACModel at (0.5, ..., 0.5)');
evaluate(HACModel, 0.5 * ones(1, getdimension(HACModel)))

disp('Prob{(U_1, ..., U_7) in the hypercube ((0.5, ..., 0.5), (0.9, ..., 0.9)]}');
prob(HACModel, 0.5 * ones(1, 7), 0.9 * ones(1, 7))

disp('The survival copula of HACModel at (0.5, ..., 0.5)');
evalsurv(HACModel, 0.5 * ones(1, 7))

%% Sample and plot data
% make the results repeatable
rng('default'); 
rng(1); % set the seed

UKnown = rnd(HACModel, 500); % do sampling
plotbimargins(UKnown);


%% Estimation
% turn the observations to pseudo-observations
U = pobs(UKnown);
% compute three estimates assuming different sets of the underlying families
fitC1219 = HACopulafit(U, getfamilies(HACModel));
fitC12 = HACopulafit(U, {'C', '12'});
fitC = HACopulafit(U, {'C'});

plot(fitC1219);
plot(fitC12);
plot(fitC);


%% Goodness-of-fit & other statistics

disp('goodness-of-fit');
[gofdSnE(fitC1219, U) ...
 gofdSnE(fitC12, U) ...
 gofdSnE(fitC, U)]

% compute an approximation of P-value 
disp('computing p value...');
tic;
estimator1 = @(U) HACopulafit(U, getfamilies(HACModel));
computepvalue(fitC1219, U, estimator1, 100)
toc

% compute the Kendall correlation matrix
K = corr(U,'type','kendall');

disp('kendall (HAC vs sample)');
[distance(fitC1219, K) ...
 distance(fitC12, K) ...
 distance(fitC, K)]

DISTANCE_TYPE = {'kendall', 'upper-tail', 'lower-tail'};
for i = 1:3
    disp([DISTANCE_TYPE{i} ' (HAC vs HAC)']);
    [distance(fitC1219, HACModel, DISTANCE_TYPE{i}) ...
     distance(fitC12, HACModel, DISTANCE_TYPE{i}) ...
     distance(fitC, HACModel, DISTANCE_TYPE{i})]
end

getdependencematrix(HACModel, 'kendall')
getdependencematrix(fitC1219, 'kendall')

getdependencematrix(HACModel, 'tails')
getdependencematrix(fitC, 'tails')

%% Auxiliaries

comparestructures(HACModel, fitC1219)
[isSameStruc, isSameFams, nSameFams] = comparestructures(fitC1219, fitC12)
tolatex(HACModel, 'cdf')
