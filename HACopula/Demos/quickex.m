% The quick example 
%
%
% Copyright 2018 Jan Gorecki

%% Installation (needed only for Octave)
% uncomment and execute from the folder 'HACopula'
% addpath(pwd);
% addpath([pwd '\' 'Demos']);
% addpath([pwd '\' 'Auxiliary']);
% addpath([pwd '\' 'Sampling']);
% savepath;
%
% For Octave 4.4 (and later versions): 
% 1) Download the 'statistics' package from 
% https://octave.sourceforge.io/statistics/index.html
% 2) Install and load the package:
% pkg install statistics-1.4.0.tar.gz 
% pkg load statistics
%
% For all versions of Octave, install similarly the symbolic package
% OctSymPy from https://github.com/cbm755/octsympy

%% Construct & plot a HAC model

% define four Archimedean generators
lam8RightRight  = {'12', tau2theta('12', 0.8)};
lam9Left        = {'19', tau2theta('19', 0.7)};
lam10Right      = {'12', tau2theta('12', 0.5)};
lam11Root       = {'C', tau2theta('C', 0.2)};

% define a 7-variate HAC
HACModel = HACopula({lam11Root, {lam9Left, 2, 5, 6}, ...
                    {lam10Right, 1, {lam8RightRight, 3, 4, 7}}});

% plot two visualizations of the HAC
plot(HACModel);
plotbipdfs(HACModel);

%% Computing probabilities involving a HAC

disp('The CDF of HACModel evaluated at (0.5, ..., 0.5)');
cdf(HACModel, 0.5 * ones(1, HACModel.Dim))

disp('Prob{(U_1, ..., U_7) in the hypercube ((0.5, ..., 0.5), (0.9, ..., 0.9)]}');
prob(HACModel, 0.5 * ones(1, HACModel.Dim), 0.9 * ones(1, HACModel.Dim))

disp('The survival copula of HACModel at (0.5, ..., 0.5)');
evalsurv(HACModel, 0.5 * ones(1, HACModel.Dim))

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
fitC     = HACopulafit(U, {'C'});

plot(fitC1219);
plot(fitC);


%% Goodness-of-fit & other statistics

disp('goodness-of-fit');
[gofdSnE(fitC1219, U) gofdSnE(fitC, U)]

% compute an approximation of p values 
disp('Computing the p value for the first estimate...');
tic;
estimatorC1219 = @(U) HACopulafit(U, getfamilies(HACModel));
computepvalue(fitC1219, U, estimatorC1219, 100)
toc
disp('Computing the p value for the second estimate...');
tic;
estimatorC = @(U) HACopulafit(U, {'C'});
computepvalue(fitC, U, estimatorC, 100)
toc


% compute the matrix of pairwise Kendall's taus
K = kendallTauMatrix(U);

% compute a distance between a HAC estimate and K
disp('kendall (HAC vs sample)');
[distance(fitC1219, K) distance(fitC, K)]

% compute a distance between a HAC estimate and the HAC model
DISTANCE_TYPE = {'kendall', 'upper-tail', 'lower-tail'};
for i = 1:3
    disp([DISTANCE_TYPE{i} ' (HAC vs HAC)']);
    [distance(fitC1219, HACModel, DISTANCE_TYPE{i}) ...
     distance(fitC, HACModel, DISTANCE_TYPE{i})]
end

% show the matrix of pairwise Kendall's tau for a HAC
getdependencematrix(HACModel, 'kendall')
getdependencematrix(fitC1219, 'kendall')

% show the matrix of upper- and lower-tail dependence coefficients a HAC
getdependencematrix(HACModel, 'tails')
getdependencematrix(fitC, 'tails')

%% Auxiliaries

% show if the structure of two HACs are the same
comparestructures(HACModel, fitC1219)

% compare HACModel to a 7-AC and get the level of match of the structures
[isSameStruc, ratioStruc] = comparestructures(HACModel, ...
                            HACopula({{'C', 0.5}, 1, 2, 3, 4, 5, 6, 7}))
% Note that the 7-AC cell model {{'C', 0.5}, 1, 2, 3, 4, 5, 6, 7} can be
% also written as [{{'C', 0.5}}, [num2cell(1:HACModel.Dim)]] , which is
% less easy to read but easier to scale.

% compare the families of two HACs
[isSameFams, ratioFams] = comparefamilies(fitC1219, fitC)

% get a LaTeX formula of a part of the HACModel's cummulative distribution
% function
tolatex(HACModel.Child{2}, 'cdf')

% access the bivariate margin of HACModel corresponding to 
% variables 1 and 3
biMargin = getbimargin(HACModel, 1, 6);

% evaluate the bivariate margin's pdf at (0.9, 0.5)
ACpdf(biMargin.Family, biMargin.Parameter, 0.9, 0.5)

