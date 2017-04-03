% The elaborated example 
%
%
% Copyright 2017 Jan Górecki

LAM8 = {'12', tau2theta('12', 0.8)};
LAM9 = {'19', tau2theta('19', 0.7)};
LAM10 = {'12', tau2theta('12', 0.5)};
LAM11 = {'C', tau2theta('C', 0.2)};
HACModel = HACopula({LAM11, {LAM9, 2, 5, 6}, {LAM10, 1, {LAM8, 3, 4, 7}}})

HACModel.Child{2}

%% Sampling a HAC elaborated
rng('default'); 
rng(1); % set the seed

UKnown = rnd(HACModel, 500);
UKnown(:, [1 2 7]) = 1 - UKnown(:, [1 2 7]);
plotbimargins(UKnown);

KNeg = corr(UKnown, 'type', 'kendall');
toFlip = findvars2flip(KNeg)

UKnown(:, toFlip) = 1 - UKnown(:, toFlip);

%% Estimating a HAC elaborated

U = pobs(UKnown);
families = getfamilies(HACModel);
[fit, fitLog] = HACopulafit(U, families, ...
                      'HACEstimator', 'pairwise',... 
                      'ThetaEstimator', 'invtau', ...
                      'ThetaEstimator2', 'invtau', ...
                      'g_1', 'average', 'g_2', @(t)mean(t), 'GOF', 'R', ...
                      'PreCollapse', true, 'Reestimator', 'Ktauavg', ...
                      'nForks', 'unknown', 'Attitude', 'optimistic', ...
                      'CheckData', 'on');

%% collapsing
% estimate binary structure
fit2Bin = HACopulafit(U, {'?'}, 'PreCollapse', false);
plot(fit2Bin);

% collapse it
K = corr(U, 'type', 'kendall');
[colHACArray, minDistArray] = collapse(fit2Bin, 'invtau', ...
                                       K, U, @(t) mean(t), ...
                                       'optimistic', 'Ktauavg', false)

% figure;
% plot(minDistArray);
% xlabel('i');
% ylabel('minDistArray(i)');

iJump = findjump(minDistArray)

fit2UnknownFams = colHACArray{iJump};
fit2 = HACopulafit(U, families, 'PreCollapsedHAC', fit2UnknownFams);
plot(fit2);
