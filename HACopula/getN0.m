function N0 = getN0(families)
%GETN0 - Get the parameter ranges of given Archimedean families
%
% N0 = getN0(families) returns a cell array N0 containing the parameter
% ranges of the input families. Note that the ranges correspond to tau >= 0
% (which is related to the SNC).
%
% References:
% [Górecki et al., 2016b] Górecki, J., Hofert, M., and Holeòa, M. (2016). On
%     structure, family and parameter estimation of hierarchical
%     Archimedean copulas. arXiv preprint arXiv:1611.09225.
% [Nelsen, 2006] Nelsen, R. (2006). An Introduction to Copulas. Springer,
%    2nd edition.
%
%
% Copyright 2017 Jan Górecki

% homogenous case
if length(families) == 1
    family = families{1};
    switch family
        case {'G', 'J', '12'}
            N0 = {{family, [1 Inf]}};
        case {'C', 'F', '19', '20'}
            N0 = {{family, [eps(0) Inf]}};
        case 'A'
            N0 = {{family, [0 1-eps(1)]}};
        case '14'
            error('HACopulafit::getN0: Due to the SNC, 14 cannot be used for homogenous HAC estimation. Use the family 14 in combination (at least) with the family C.');
        case '?'
            N0 = {{family, [-1 1]}};  % an arbitrary interval that corresponds to the maximal range of Kendall's tau          
        otherwise
            error(['HACopulafit::getN0: Unsupported family - ', family]);
    end
    return;
end
   
% heterogenous case

%F24 = {'C', '12', '14', '19', '20'};
%F1234 = {'A', 'C', '19', '20'};

if      ((sum(strcmp('A',families)) == 1) && (sum(strcmp('C',families)) == 1)) || ...
        ((sum(strcmp('A',families)) == 1) && (sum(strcmp('20',families)) == 1)) || ...
        ((sum(strcmp('A',families)) == 1) && (sum(strcmp('19',families)) == 1))
    % the case when a pair from \mathcal{L}_3 ([Górecki et al., 2016b]) is provided or
    % only the pair from \mathcal{L}_1 ([Górecki et al., 2016b]) is provided
    
    N0 = { {'A', [0 1-eps(1)]}, {'C', [1 Inf]}, {'19', [eps(0) Inf]}, {'20', [1 Inf]} };
    
elseif  ((sum(strcmp('C',families)) == 1) && (sum(strcmp('12',families)) == 1)) || ...
        ((sum(strcmp('C',families)) == 1) && (sum(strcmp('19',families)) == 1)) || ...
        ((sum(strcmp('C',families)) == 1) && (sum(strcmp('14',families)) == 1)) || ...
        ((sum(strcmp('C',families)) == 1) && (sum(strcmp('20',families)) == 1))
    % the case when a pair from \mathcal{L}_2 or \mathcal{L}_4 ([Górecki et al., 2016b]) is provided
    
    N0 = { {'C', [eps(0) Inf]}, {'12', [1 Inf]}, {'14', [1 Inf]},...
        {'19', [eps(0) Inf]}, {'20', [eps(0) Inf]} };
    
else
    error(['HACopulafit::getN0: The family set {' strjoin(families, ', ') '} is not supported for the input ''families''.' ...
           ' The input parameter ''families'' should be either '...
            ' 1) a subset of  {A, C, 19, 20} provided it contains the family A, or 2) a subset of {C, 12, 14, 19, 20}' ...
            ' provided it contains the family C. Sets of families that do not contain neither A nor C are not supported,' ...
            ' as there is no guarantee following from the SNC that a proper copula results.']);
end



N0families = [N0{:}];
N0families = N0families(1:2:end);

% check if families is a subset of  N0 families
for i = 1:length(families)
    if (sum(strcmp(families(i), N0families)) == 0)
            error(['HACopulafit::getN0: The family set {' strjoin(families, ', ') '} is not supported for the input ''families''.' ...
           ' The input parameter ''families'' should be either '...
            ' 1) a subset of  {A, C, 19, 20} provided it contains the family A, or 2) a subset of {C, 12, 14, 19, 20}' ...
            ' provided it contains the family C. Sets of families that do not contain neither A nor C are not supported,' ...
            ' as there is no guarantee following from the SNC that a proper copula results.']);
    end
end

keepInN0 = zeros(1, length(N0));
for i = 1:length(N0)
    keepInN0(i) = (sum(strcmp(N0families(i), families)) == 1);
end
% remove what is not needed
N0 = N0(keepInN0 == 1);


