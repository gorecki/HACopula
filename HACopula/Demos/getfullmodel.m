function model = getfullmodel(nLevels, dAC, family, tauRoot, tauDiff)
%GETFULLMODEL - construct a type of a HAC resembling a fully-nested AC
%
% Returns a HAC model with all generators from FAMILY
% with nesting levels equal to NLEVELS and one DAC-dimensional AC at
% each level, i.e., d = nLevels * dAC - (nLevels-1), such that tau =
% TAUROOT for the root and tau(child)-tau(parent) = TAUDIFF.
%
% Example:
% getfullmodel(11, 10, 'C', 0.1, 0.08) returns a 100-HAC from the Clayton
% family with tau = 0.1 at the root and tau(child)-tau(parent) = 0.08.
%
%
% Copyright 2018 Jan Gorecki

%% HACopula cell structure construction
% the lowest level (without generator)
HACCellModel = num2cell((nLevels-1)*(dAC-1)+1:nLevels*(dAC-1)+1);
HACCellModel = [{{family, tau2theta(family ,(nLevels-1) * tauDiff + tauRoot) }} HACCellModel];

% construct HACopula cell structure
for iLevel = nLevels-1:-1:1
    ACCellModel = num2cell((iLevel-1)*(dAC-1)+1:iLevel*(dAC-1));
    % add generator
    ACCellModel = [{{family, tau2theta(family,(iLevel-1) * tauDiff + tauRoot)}} ACCellModel]; 
    HACCellModel = [ACCellModel {HACCellModel}];
end

model = HACopula(HACCellModel);

end