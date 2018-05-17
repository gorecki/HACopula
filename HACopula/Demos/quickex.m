% The quick example 
%
%
% Copyright 2017 Jan Górecki

%% Installation (needed only for Octave)
% uncomment and execute from the folder 'HACopula'
% addpath(pwd);
% addpath([pwd '\' 'Demos']);
% addpath([pwd '\' 'Auxiliary']);
% addpath([pwd '\' 'Sampling']);
% savepath;

%% Construct & plot a HAC model

% define four Archimedean generators
LAM8 = {'12', tau2theta('12', 0.8)}
LAM9 = {'19', tau2theta('19', 0.7)}
LAM10 = {'12', tau2theta('12', 0.5)}
LAM11 = {'C', tau2theta('C', 0.2)}

% define a 7-variate HAC
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

% make the results repeatable (as the random number generator of Octave
% differs from the one in MATLAB (see
% https://octave.sourceforge.io/octave/function/rand.html), if we are in
% Octave, load the data (UKnown) sampled in MATLAB in order to get the
% further results same for both Octave and MATLAB)
if isoctave % are we in Octave?
    load('Demos/quickex_data.mat'); 
else % we are in MATLAB
    % set the seed
    rng('default'); 
    rng(1);
    % sample 500 random vectors from HACModel
    UKnown = rnd(HACModel, 500); 
end

plotbimargins(UKnown);


%% Estimation

% turn the observations (the sample from HACModel) to pseudo-observations
U = pobs(UKnown);
% compute three HAC estimates assuming different sets of the underlying families
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

% compute an approximation of p value 
disp('computing p value...');
tic;
estimator1 = @(U) HACopulafit(U, getfamilies(HACModel));
computepvalue(fitC1219, U, estimator1, 100)
toc

% compute the Kendall correlation matrix
K = kendallTauMatrix(U);

% compute a distance between a HAC estimate and K
disp('kendall (HAC vs sample)');
[distance(fitC1219, K) ...
 distance(fitC12, K) ...
 distance(fitC, K)]

% compute a distance between a HAC estimate and the HAC model
DISTANCE_TYPE = {'kendall', 'upper-tail', 'lower-tail'};
for i = 1:3
    disp([DISTANCE_TYPE{i} ' (HAC vs HAC)']);
    [distance(fitC1219, HACModel, DISTANCE_TYPE{i}) ...
     distance(fitC12, HACModel, DISTANCE_TYPE{i}) ...
     distance(fitC, HACModel, DISTANCE_TYPE{i})]
end

% show the matrix of Kendall's tau for a HAC
getdependencematrix(HACModel, 'kendall')
getdependencematrix(fitC1219, 'kendall')

% show the matrix of upper- and lower-tail dependence coefficients a HAC
getdependencematrix(HACModel, 'tails')
getdependencematrix(fitC, 'tails')

%% Auxiliaries

% show if the structure of two HACs are the same
comparestructures(HACModel, fitC1219)
% also compare the families of two HACs
[isSameStruc, isSameFams, nSameFams] = comparestructures(fitC1219, fitC12)
% get a LaTeX formula of a part of the HACModel's cummulative distribution
% function
tolatex(HACModel.Child{2}, 'cdf')
