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
if isoctave
    load('Demos/quickex_data.mat'); 
else % MATLAB    
    % set the seed
    rng('default'); 
    rng(1);
    % do sampling
    UKnown = rnd(HACModel, 500); 
end

% introduce some negative correlation by flipping three variables (columns)
UKnown(:, [1 2 7]) = 1 - UKnown(:, [1 2 7]);
plotbimargins(UKnown);

KNeg = kendallTauMatrix(UKnown);

% find variables that can reduce the negative correlation if they are
% flipped
toFlip = findvars2flip(KNeg)
% flip the variables identified in the previous step
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

%% Collapsing a binary structured HAC

% fit a binary structured HAC assuming unknown ('?') families
fit2Bin = HACopulafit(U, {'?'}, 'PreCollapse', false);
plot(fit2Bin);

K = kendallTauMatrix(U);
% collapse the binary structured HAC
[colHACArray, minDistArray] = collapse(fit2Bin, 'invtau', ...
                                       K, U, @(t) mean(t), ...
                                       'optimistic', 'Ktauavg', false)

% show the minimal distances from the collapsing process
figure;
plot(minDistArray);
xlabel('i');
ylabel('minDistArray(i)');

% find the first substantial jump
iJump = findjump(minDistArray)

% take the HAC corresponding to this jump
fit2UnknownFams = colHACArray{iJump};
% estimate the families and parameters assuming the structure of fit2UnknownFams
fit2 = HACopulafit(U, families, 'PreCollapsedHAC', fit2UnknownFams);
plot(fit2);
