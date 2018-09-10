function [collapsedHACArray, minDistanceArray] = collapse(obj, ...
    thetaEstimatorPairwise, K, U, g,...
    attitude, reEstimateType, considerFamilies)
%COLLAPSE - return an array of HACs simplified by collapsing
% two forks into one repeatedly
%
% Usage:
% [collapsedHACArray, minDistanceArray] = collapse(obj, thetaEstimatorPairwise, K, U, g, attitude, reEstimateType, considerFamilies)
%
% Purpose: The function repeatedly collapses two forks into one
% until there is only one fork remaining (i.e., until the
% original HAC collapses to an AC). In each collapsing step,
% two parent-child forks with the lowest absolute difference between their
% taus are found and then collapsed into one fork. This
% collapsed fork inherits the family of the original forks, and
% its parameter is re-estimated using a provided approach to
% re-estimation (Ktauavg or taumin). Also, in order to the
% SNC be satified for this fork, the re-estimated parameter
% is trimed (for the optimistic attitude) to an appropriate
% parameter range.  Then, the Tau property is re-computed
% according to the new parameter value and also TauOrdering of
% the whole collapsed HACopula object is re-calculated.
%
% Inputs:
% obj                    - A HACopula object to collapse.
% thetaEstimatorPairwise - An estimator for the paramter of a
%                          generator ('invtau', 'invtau2',
%                          'mle' or 'mle2').
% K                      - The Kendall correlation matrix of U.
% U                      - The observations that were used for
%                          estimation of obj.
% g                      - A [0, 1]-aggregation function from
%                          {'average', 'min', 'max'}.
% attitude               - 'optimistic' or 'pessimistic'.
% reEstimateType         - 'Ktauavg' - the [Gorecki et al.,
%                          2016b] approach,
%                          or 'taumin' - the [Uyttendaele, et
%                          al., 2016] approach.
% considerFamilies       - true or false.
%                          If true, only
%                          parent-child pairs from the same
%                          family are considered for
%                          collapsing. This also implies that
%                          the collapsing process stops if all
%                          parent-child families in the actual
%                          HAC are different, i.e., the
%                          collapsing to an AC is not reached.
%                          If false (default), the families are
%                          ignored. Note that this option
%                          enables to successfully use the
%                          function findjump, whereas for the
%                          option considerFamilies = true,
%                          findjump generally fails to return
%                          a proper number of forks in the true
%                          HAC.
%
% Outputs:
% collapsedHACArray      - A cell array of HACopula objects
%                          simplified from obj by collapsing.
%                          Note that collapsedHACArray{1}
%                          contains the original obj,
%                          collapsedHACArray{2} contains a HAC
%                          with one fork less then in obj, etc.
% minDistanceArray       - An array of minimal distances. Note
%                          that minDistanceArray(1) = 0.
%                          minDistanceArray(2) = |tau1 - tau2|,
%                          where tau1 and tau2 correspond to
%                          the taus of the two forks in
%                          collapsedHACArray{1} that are
%                          collapsed into one forks in
%                          collapsedHACArray{2}, etc.
%
% NOTE: The re-estimation of the parameter of the collapsed
% node is performed only by using the pairwise estimation. This
% is due to the fact that MLE needs densities of ACs for d >
% 2, which are however uknown for the concerned AC families
% (except (A, C, F, G, J)). Also, because of such an approach
% to the re-estimation, U and K that were used for estimation
% of the HACopula input obj are needed.
%
%
% Copyright 2018 Jan Gorecki

if ~(strcmp(reEstimateType, 'Ktauavg') || strcmp(reEstimateType, 'taumin'))
    error('HACopula:collapse', 'HACopula::collapse: the parameter reEstimateType must be  either ''Ktauavg'' or ''taumin''.')
end

nForks = size(obj.Forks, 2);

if nForks == 1
    error('HACopula:collapse', 'HACopula::collapse: Attempting to collapse HAC that is AC.')
end

collapsedHACArray = cell(1, nForks);
collapsedHACArray{1} = obj; % return also the original HAC

minDistanceArray = zeros(1, nForks); % store the minimum distances for each collapsed HAC
% minDistanceArray(1) = 0;  % the minimal distance of the original HAC is kept to be 0

families = getfamilies(obj);

for k = 2:nForks
    
    cpObj = copy(collapsedHACArray{k-1});
    
    minDistance = Inf;
    
    % find the parent-child pair (parentToCollapse, childToCollapse) to collapse
    for j = 1:size(cpObj.Forks, 2)
        parent = cpObj.Forks{j};
        for i = 1:length(parent.Child)
            child = parent.Child{i};
            if isa(child, 'HACopula') && (strcmp(parent.Family, child.Family) || ~considerFamilies)
                distance = child.Tau - parent.Tau; % compute the distance
                if distance < minDistance
                    parentToJoin = parent;
                    childToJoin = child;
                    iChild = i;
                    minDistance = distance;
                end
            end
        end
    end
    
    if isinf(minDistance)
        % there are no more forks to collapse due to different
        % families in forks -> stop collapsing
        collapsedHACArray = collapsedHACArray(1:k-1);
        minDistanceArray = minDistanceArray(1:k-1);
        return
    end
    
    % join the pair in the HAC
    nChildren = size(parentToJoin.Child, 2);
    % assing all (except the first) the children of the childToJoin to parentToJoin
    for j = 2:size(childToJoin.Child, 2);
        if isa(childToJoin.Child{j}, 'HACopula')
            decreaselevel(childToJoin.Child{j});
            childToJoin.Child{j}.Parent = parentToJoin;
        end
        parentToJoin.Child{nChildren + j - 1} = childToJoin.Child{j};
    end
    % replace the childToJoin by its first child
    if isa(childToJoin.Child{1}, 'HACopula')
        decreaselevel(childToJoin.Child{1});
        childToJoin.Child{1}.Parent = parentToJoin;
    end
    parentToJoin.Child{iChild} = childToJoin.Child{1};
    
    % remove childToJoin from whole HAC (from forks)
    % first from the root
    for j = 1:size(cpObj.Forks,2)
        if cpObj.Forks{j}.TauOrdering == childToJoin.TauOrdering
            remaining = setdiff(1:size(cpObj.Forks,2), j);
            cpObj.Forks = cpObj.Forks(remaining);
            break
        end
    end
    % then from the rest of the forks
    for i = 2:size(cpObj.Forks, 2)
        fork = cpObj.Forks{i};
        for j = 1:size(fork.Forks,2)
            if fork.Forks{j}.TauOrdering == childToJoin.TauOrdering
                remaining = setdiff(1:size(fork.Forks,2), j);
                fork.Forks = fork.Forks(remaining);
                break
            end
        end
    end
    
    % finally, delete the childToJoin fork
    delete(childToJoin);
    
    % now, re-estimate the parameter of parentToJoin
    if strcmp(reEstimateType, 'Ktauavg')
        % use the re-estimation proposed in [Gorecki et al., 2017]
        % store the descendant leaves
        nChildren = size(parentToJoin.Child, 2);
        descLeaves = cell(1, nChildren);
        for i = 1:nChildren
            if isa(parentToJoin.Child{i}, 'HACopula')
                descLeaves{i} = parentToJoin.Child{i}.Leaves;
            else
                descLeaves{i} = parentToJoin.Child{i};
            end
        end
        % trim the parameter in order to satisfy the SNC
        % (1) first do trim to satisfy it with all children
        N = getN0(families);
        for i = 1:nChildren
            if isa(parentToJoin.Child{i}, 'HACopula')
                N = intersectNs(N, computeN2(parentToJoin.Child{i}.Family, parentToJoin.Child{i}.Parameter, families));
            end
        end
        % get the families that remainded in N
        Nfamilies = [N{:}];
        Nfamilies = Nfamilies(1:2:end);
        % find the index of the parentToJoin family
        iFamily = strcmp(Nfamilies, parentToJoin.Family);
        % get its SNC range
        r = N{iFamily}{2};
        % adjust the parentToJoin's parameter (note that the family of the parent and of the child is the same)
        parentToJoin.Parameter = fitparameter(descLeaves, parentToJoin.Family, 0, thetaEstimatorPairwise, '', K, U, g, [], [], parentToJoin.Tau, r, attitude); % NOTE: uses only the pairwise estimator
        % trim the parameter
        parentToJoin.Parameter = min(max(r(1),parentToJoin.Parameter),r(2));
        % (2) then do trim to satisfy SNC with the parent of parentToJoin
        if parentToJoin.Level > 1
            if strcmp(parentToJoin.Parent.Family, parentToJoin.Family)
                % homogenous case
                parentToJoin.Parameter = max(parentToJoin.Parameter, parentToJoin.Parent.Parameter);
            else
                % heterogenous case
                switch [parentToJoin.Parent.Family parentToJoin.Family]
                    %case 'A19'                                          % the L_1 class
                    % no trim necessary
                    %case {'C12', 'C19'}                                 % the L_2 class
                    % no trim necessary
                    case {'AC', 'A20'}                                  % the L_3 class
                        parentToJoin.Parameter = max(parentToJoin.Parameter, 1); % need to be at least 1
                        %case 'C14'
                        % two families of 14 are never going to
                        % collapse as the SNC is not known
                        % parentToJoin.Parameter = min(parentToJoin.Parameter, 1/parentToJoin.Parent.Parameter);
                    case 'C20'                                          % the L_4 class
                        parentToJoin.Parameter = max(parentToJoin.Parameter, parentToJoin.Parent.Parameter);
                end
            end
        end
    else % if strcmp(reEstimateType, 'taumin')
        % use the re-estimation proposed in [Uyttendaele, et al., 2016]
        % i.e. use the mimimum of the tau equivalents of the parameters of the children, i.e., keep the original parent parameter
        % parentToJoin.Parameter = parentToJoin.Parameter; % i.e, do nothing
    end
    % re-compute tau
    parentToJoin.Tau = theta2tau(parentToJoin.Family, parentToJoin.Parameter);
    
    % NOTE: necessary only for the plot function (comment to
    % make it faster)
    addtauordering(cpObj);
    
    collapsedHACArray{k} = cpObj;
    
    % store the minimal distance
    minDistanceArray(k) = minDistance;
end
end
