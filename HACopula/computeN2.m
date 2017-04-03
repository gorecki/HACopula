function N2 = computeN2(family, theta, families)
%COMPUTEN2 - Computing admissible parameter ranges for the parent of its
% child generator.
%
% N2 = computeN2(family, theta, families) returns a cell array containing admissible parameter ranges for the 
% generator from family with the parameter theta also considering all
% possible families. This function implements the mapping
% \mathcal{N}^2_{\mathcal{F}}(a, theta) introduced in [Górecki et al., 2016b], where
% \mathcal{F} there corresponds to the input families here and the argument
% a there corresponds to the parameter family here.
%
% Example:
% N2 = computeN2('20', 1.5, {'A', 'C', '19', '20'}) returns the following three
% parameter ranges: [0, 1) for the family A, (0, 1.5] for the family C and
% (0, 1.5] for the family 20. Observe that any parent generator from these
% three families with the parameter from the corresponding parameter range would
% satisfy the sufficient nesting condition with the child generator from the
% family 20 with the parameter 1.5. In other words, computeN2 returns
% parameter ranges, from which the parameter of parent generator should be
% selected in the estimation process in order to satisfy the sufficient
% nesting condition. Note that the input families defines all families from
% which the parent generator could be selected (but before considering the
% sufficient nesting condition) - this is why the family 19 is not in the output.
%
% References:
% [Górecki et al., 2016b] Górecki, J., Hofert, M., and Holeòa, M. (2016). On
%     structure, family and parameter estimation of hierarchical
%     Archimedean copulas. arXiv preprint arXiv:1611.09225.
%
%
% Copyright 2017 Jan Górecki

% homogenous case
if length(families) == 1
    family = families{1};
    switch family
        case {'G', 'J', '12', '14'}
            N2 = {{family, [1 theta]}};
        case {'C', 'F', '19', '20'}
            N2 = {{family, [eps(0) theta]}};
        case 'A'
            N2 = {{family, [0 min(1-eps(1), theta)]}};
        case '?'
            N2 = {{family, [-1 1]}};
        otherwise
            error(['HACopulafit::computeN2: Unsupported family - ', family]);
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
    
    switch family
        case 'A'
            N2 = {{'A', [0 theta]}};
        case 'C'
            N2 = {{'A', [0 1-eps(1)]}, {'C', [eps(0) theta]}};
        case '19'
            N2 = {{'A', [0 1-eps(1)]}, {'C', [eps(0) 1]}, {'19', [eps(0) theta]}};
        case '20'
            N2 = {{'A', [0 1-eps(1)]}, {'C', [eps(0) theta]}, {'20', [eps(0) theta]}};
        otherwise
            error(['HACopulafit::computeN2: unsupported family ' family '.']);
    end
    
    
elseif  ((sum(strcmp('C',families)) == 1) && (sum(strcmp('12',families)) == 1)) || ...
        ((sum(strcmp('C',families)) == 1) && (sum(strcmp('19',families)) == 1)) || ...
        ((sum(strcmp('C',families)) == 1) && (sum(strcmp('14',families)) == 1)) || ...
        ((sum(strcmp('C',families)) == 1) && (sum(strcmp('20',families)) == 1))
    % the case when a pair from \mathcal{L}_2 or \mathcal{L}_4 ([Górecki et al., 2016b]) is provided
    
    % compute the output of \mathcal{N}^2_{\mathcal{F}_{24}}
    switch family
        case 'C'
            N2 = {{'C', [eps(0) theta]}};
        case '12'
            N2 = {{'C', [eps(0) 1]}, {'12', [1 theta]}};
        case '14'
            N2 = {{'C', [eps(0) 1/theta]}};
        case '19'
            N2 = {{'C', [eps(0) 1]}, {'19', [eps(0) theta]}};
        case '20'
            N2 = {{'C', [eps(0) theta]}, {'20', [eps(0) theta]}};
        otherwise
            error(['HACopulafit::computeN2: unsupported family ' family '.']);
    end
else
    error('HACopulafit::computeN2: Unsupported family combination in the parameter ''families''');
end

% remove the families that are not in the parameter 'families'
N2families = [N2{:}];
N2families = N2families(1:2:end);
keepInN2 = zeros(1, length(N2));
for i = 1:length(N2)
    keepInN2(i) = (sum(strcmp(N2families(i), families)) == 1);
end
% remove what is not needed
N2 = N2(keepInN2 == 1);

