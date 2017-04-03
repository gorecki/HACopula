classdef HACopula < matlab.mixin.Copyable
    %HACopula - The class of hierarchical Archimedean copula (HAC) models
    %
    % As a HAC can be viewed as a set of generators (representing
    % Archimedean copulas) that are recursively nested into each other, an
    % object from the HACopula class is also defined recursively. This
    % means, instantiating a HACopula object, it contains a description of
    % the generator at the root of the HAC structure, which is given by the
    % family and the parameter of the generator, and a cell array of
    % HACopula objects or integers, which represent the children (another
    % HACs or leaves) of a generator in the HAC structure. Hence the
    % three essential properties of a HACopula object are Family,
    % Parameter and Child. The rest of the properties are closely related
    % to (and hence useful for) HAC models (Level, Tau, TauOrdering) or are
    % used for increasing computational efficiency (Leaves, Parent, Root,
    % Forks).
    %
    % Note that HACopula is a handle class, i.e., modifications of a
    % HACopula object passed as an input to a function that modifies this
    % input are propagated to the input object. If a method is modifying
    % its HACopula input object, this is stressed out in the help of this
    % function. To prevent a HACopula object from modification in such
    % functions, use the copy method and pass the copied object to these
    % functions.
    %
    % Also note that we will sometimes use the term "fork" instead of
    % "generator", as there is a clear relationship between these two,
    % e.g., see [Górecki et al., 2016b].
    %
    %
    % References:
    % [Genest et al., 2009] Genest, C., Rémillard, B., and Beaudoin, D.
    %    (2009). Goodness-of-fit tests for copulas: A review and a power
    %    study. Insurance: Mathematics and Economics, 44(2):199-213.
    % [Górecki et al., 2014] Górecki, J., Hofert, M., and Holeòa, M. (2014). On
    %     the consistency of an estimator for hierarchical Archimedean copulas.
    %     In 32nd International Conference on Mathematical Methods in
    %     Economics, pages 239-244.
    % [Górecki et al., 2016b] Górecki, J., Hofert, M., and Holeòa, M. (2016). On
    %     structure, family and parameter estimation of hierarchical
    %     Archimedean copulas. arXiv preprint arXiv:1611.09225.
    % [Nelsen, 2006] Nelsen, R. (2006). An Introduction to Copulas. Springer,
    %     2nd edition.
    % [Uyttendaele, et al., 2016] Uyttendaele, et al., On the estimation of
    %     nested Archimedean copulas: A theoretical and an experimental
    %     comparison, Tech. rep., UCL (2016).
    %
    %
    % Copyright 2017 Jan Górecki and Martin Holeòa
    
    properties
        Family          % A family label from {'A', 'C', 'F', 'G', 'J', '12', '14', '19', '20', '?'}, where
                        % A = Ali-Mikhal-Haq, C = Clayton,
                        % F = Frank,
                        % G = Gumbel,
                        % J = Joe, and
                        % 12, 14, 19 and 20 - see page 116 in [Nelsen,
                        % 2006].
                        % The parametric forms of the families listed above
                        % can be found in the function getgenerator. The family '?' is an
                        % auxiliary family that is used in cases, when one
                        % wants to estimate (or also to collapse)
                        % the structure of a HAC without any assumptions on
                        % the generators of the underlying HAC. This family
                        % thus cannot be used for purposes of parameter
                        % estimation. Also note that, to allow this
                        % family to be a valid family of HACopula, its
                        % parameter is arbitrarily defined to be equal to
                        % the Kendall's tau.
        Parameter       % A real value >= 0.
        Tau             % The value of Kendall's tau corresponding to Parameter 
                        % and Family.
        TauOrdering     % The order of the generator in the tau ordering 
                        % defined by the function addtauordering. Serves as an identifier of of the fork.
        Level           % The nesting level of a generator (Level = 1 for the root).
        Leaves          % A cell array containing all descendant leaves of the actual fork.
        % HACopula handles
        Child           % A cell array of child (HACopula) objects.
        Parent          % [], if Level = 1, otherwise it contains the (HACopula) parent of this object.
        Root            % The (HACopula) root of this  object.
        Forks           % A cell array with HACopula objects containing all 
                        % descendant forks of this fork including this
                        % fork.
    end
    
        methods(Access = protected)
            

            function cpObj = copyElement(obj)
                % returns a copy of a HACopula object
                % Override copyElement method:
                
                % Make a shallow copy of all properties
                cpObj = copyElement@matlab.mixin.Copyable(obj);
                % redirect to the copied forks
                cpObj.Forks = {cpObj};                
                % Make a deep copy 
                for i = 1:size(obj.Child,2)
                    if isobject(obj.Child{i})
                        cpObj.Child{i} = copy(obj.Child{i});
                        % redirect to the copied parent
                        cpObj.Child{i}.Parent = cpObj;
                        % redirect to the copied forks
                        cpObj.Forks = [cpObj.Forks cpObj.Child{i}.Forks];
                    end
                end
                % redirect to the copied root
                if cpObj.Level == 1
                    for i = 1:size(cpObj.Forks, 2)
                        cpObj.Forks{i}.Root = cpObj;
                    end
                end
            end
        end
    
    
    methods
        
        function obj = HACopula(model)
            % The constructor of the class.
            % Accepts a cell nested structure or another HACopula object.
            % Checks if leaves are {1, ..., d}.
            % Checks for the SNC.
            % Adds an ordering of the forks according to tau.
            
            if nargin == 0
                % return an empty HAC
                return;
            end
            
            narginchk(1,3);
            
            if iscell(model)
                cell2tree(obj, model);
            elseif isa(model,'HACopula') % is HACopula?
                obj = copy(model); % make a copy
            else
                error('HACopula: Unrecognized HAC model.');
            end
            
            % do basic checks
            checkleaves(obj);
            checksnc(obj);
            % add tau ordering
            addtauordering(obj);
        end
        
        
        function addtauordering(obj)
            % Adds an ordering of the forks according to Kendall's tau, which is
            % stored in the TauOrdering property. The forks are
            % ordered according to corresponding taus. The root, which
            % has always the minimal tau value, gets the position (#forks + #leaves) in this order,
            % the fork with the closest higher value of tau gets the position (#forks +
            % #leaves - 1), ..., and the fork with the maximal tau value
            % gets the position (#forks + 1). If there are two forks with the
            % same tau value, the ordering involves the lexicographical
            % ordering of their descendant leaves.
            %
            % NOTE: This function modifies the input obj.
            
            d = getdimension(obj);
            
            % get leaves of each of the forks and extend them
            nForks = size(obj.Forks,2);
            leavesExtend = zeros(nForks, d);
            tauOfForks = zeros(1, nForks);
            for i = 1:nForks
                % extend leaves for sortrow below
                leavesExtend(i, 1:d) =  [ones(1, d - size(obj.Forks{i}.Leaves,2))*(-Inf) sort(obj.Forks{i}.Leaves)];
                tauOfForks(i) = obj.Forks{i}.Tau;
            end
            
            uniqueTauOfForks = unique(tauOfForks);
            tauOrdering = d + 1;        % start from d + 1
            forkNumbers = 1:nForks;
            
            leavesTauOrdering = zeros(1, nForks);
            for iUnFork = size(uniqueTauOfForks, 2):-1:1    % start from the highest tau
                sameTaus = (tauOfForks == uniqueTauOfForks(iUnFork));
                [~,lexOrdering] = sortrows(leavesExtend(sameTaus,:));
                [~, lexOrdering] = sort(lexOrdering);           % get indices of the ordering
                
                leavesTauOrdering(forkNumbers(sameTaus)) = lexOrdering + tauOrdering - 1;   % associate each set of leaves with a tau ordering number
                tauOrdering = tauOrdering + sum(sameTaus);
            end
            
            % add tau ordering
            for i = 1:nForks
                obj.Forks{i}.TauOrdering = leavesTauOrdering(i);
            end
        end
        
        
        
        function adjacencyMatrix = adjacencymat(obj)
            % ADJACENCYMAT - a d x d matrix showing 1 if two nodes are connected in
            % HAC structure and 0, otherwise.
            %
            % Note that the forks are labeled according to TauOrdering.
            %
            % Example:
            %
            % myHAC = HACopula({{'A', .5}, 1, {{'A', 0.9}, 2, {{'C', 2}, 3, 4}}});
            % adjacencymat(myHAC)
            % 
            % ans =
            % 
            %      0     0     0     0     0     0     1
            %      0     0     0     0     0     1     0
            %      0     0     0     0     1     0     0
            %      0     0     0     0     1     0     0
            %      0     0     1     1     0     1     0
            %      0     1     0     0     1     0     1
            %      1     0     0     0     0     1     0
            
            
            orderedNodes = obj.Forks;
            nNodes = size(orderedNodes, 2) + size(obj.Leaves,2);
            adjacencyMatrix = zeros(nNodes, nNodes);
            for node = 1:length(orderedNodes)
                thisNode = orderedNodes{node};
                childArray = zeros(1, size(thisNode.Child, 2));
                for j = 1:length(thisNode.Child)
                    if isa(thisNode.Child{j}, 'HACopula')
                        childArray(j) = thisNode.Child{j}.TauOrdering;
                    else
                        childArray(j) = thisNode.Child{j};
                    end
                end
                
                adjacencyMatrix(thisNode.TauOrdering, childArray) = ones(1, size(childArray,2));
            end
            % make it symetric
            adjacencyMatrix = adjacencyMatrix + adjacencyMatrix';
        end
        
        
        
        function [bootstrapSamples,theirSn] = bootstrapSn(obj, estimator, bootstrapSize, sampleSize)
            % BOOTSTRAPSN - a bootstrap for a copula estimate
            %
            % Inputs:
            % obj           - A HACopula object.
            % estimator     - A HAC estimator with one input argument 
            %                 corresponding to pseudoobservations, e.g., an
            %                 anonymous funtion @(U)f(U).
            % bootstrapSize - A bootstrap size (default 1000).
            % sampleSize    - A sample size (default 100).
            %
            % For a given HAC, performs the bootstrap method proposed in
            % [Genest et al., 2009], returns the bootstrap samples and
            % their values of the rank-based Cramér-von Mises statistics
            % (denoted Sn in the paper). 
            
            narginchk(2,4);
            
            if nargin < 3
                bootstrapSize = 1000;
                sampleSize = 100;
            elseif nargin < 4
                sampleSize = 100;
            end
            
            % start bootstrapping process
            
            bootstrapSamples = cell(1, bootstrapSize);
            theirSn = zeros(1, bootstrapSize);
            
            for sample = 1:bootstrapSize
%                 if mod(sample,10) == 1
%                     disp(sample);
%                 end
                % sample
                UKnown = rnd(obj, sampleSize);
                % turn to pseudo-observations
                U = pobs(UKnown);
                bootstrapSamples{sample} = U;
                % compute GoF for the estimate on the sample
                theirSn(sample) = gofdSnE(estimator(U), U);
            end
        end
        
        
        
        
        function cell2tree(obj, HACCellModel, varargin)
            %CELL2TREE - convert a nested cell structure to a HACopula
            %object
            % 
            % Given a specified cell nested structure, the input obj is
            % modified in a way that it is a HACopula object describing the
            % HAC given by the specified cell nested structure.
            %
            % Defining the nested structure, at each level, first specify
            % the generator {family, parameter}, e.g., {'C', 1.5}, and then
            % its children, e.g., {{'C', 1.5}, 1, 2} corresponds to the
            % 2-AC Clayton with theta = 1.5.  Note that each child could be
            % another nested cell structure.
            %
            % Example (3 nested ACs):
            % myHAC = HACopula(); % the empty HAC model
            % cellModel = {{'A', 0.1},{{'A', 0.5}, 3, 4, 5},{{'C', 1.5}, 1, 2}};
            % cell2tree(myHAC, cellModel);
            %
            % myHAC contains the HAC model corresponding to cellModel, use
            % plot(myHAC) to see it graphically
            %
            % Example 2 (9-HAC):
            % cellModel = {{'A', 0.1},{{'A', 0.5}, {{'19', 2}, 3, 6, 7}, 4, 5}, ...
            % {{'C', 1.5}, 1, {{'20', 1.25}, 2, 8 ,9}}}
            %
            % NOTE: This function modifies the input obj.
            
            narginchk(2, 3)
            
            if size(varargin, 2) == 0
                level = 1;
                obj.Root = obj;
            else
                level = varargin{1};
            end
            
            nCells = size(HACCellModel, 2);
            
            obj.Level = level;
            
            if nCells > 2
                % it is a fork
                for i = 1:nCells
                    isFound = false;
                    % check cell structure
                    if (size(HACCellModel{i},2) == 2) && (ischar(HACCellModel{i}{1})) && ...
                            (isnumeric(HACCellModel{i}{2}))
                        % check if the first cell is a generator
                        if isgenerator(HACCellModel{i}{1}, HACCellModel{i}{2})
                            isFound = true;
                            iGenerator = i;
                            break;
                        else
                            error(['HACopula.cell2tree: The generator (' HACCellModel{i}{1} ', ' num2str(HACCellModel{i}{2}) ') is not supported.']);
                        end
                    end
                end
                if ~isFound
                    error('HACopula.cell2tree: wrong cell HAC structure.');
                end
                obj.Family = HACCellModel{iGenerator}{1};
                obj.Parameter = HACCellModel{iGenerator}{2};
                obj.Tau = theta2tau(obj.Family, obj.Parameter);
                obj.Forks = {obj}; % add this fork to Forks
                
                nChildren = nCells - 1;
                for iChild = 1:nChildren
                    if size(HACCellModel{iChild+1},2) == 1
                        % it is a leaf
                        obj.Child{iChild} = HACCellModel{iChild+1};
                        obj.Leaves = [obj.Leaves HACCellModel{iChild+1}];
                    else
                        obj.Child{iChild} = HACopula();
                        cell2tree(obj.Child{iChild}, HACCellModel{iChild+1}, level + 1);
                        obj.Child{iChild}.Parent = obj;
                        obj.Child{iChild}.Root = obj.Root;
                        obj.Leaves = [obj.Leaves obj.Child{iChild}.Leaves];
                        obj.Forks = [obj.Forks obj.Child{iChild}.Forks];
                    end
                end
            else
                error('HACopula.cell2tree: Wrong cell HAC structure. Invalid number of subcells.')
            end
        end
        
        
        
        function checkleaves(obj)
            %CHECKLEAVES - check if the leaves of a HACopula object consitute {1, ..., d}.
            
            %get leaves of the model
            leaves = obj.Leaves;
            
            % check if the leaves are integers
            if sum(floor(leaves) == leaves) ~= size(leaves,2)
                error('HACopula.checkleafs: The leaves of the HAC are not integer numbers.');
            end
            
            uniqueLeaves = unique(leaves);
            if size(uniqueLeaves,2) ~= size(leaves,2)
                error('HACopula.checkleafs: The leaves of the HAC are not unique, i.e., there are two or more of the same leaves.');
            end
            
            % the dimension of the HAC is the number of its leaves
            d = size(leaves, 2);
            
            % the leaves must be 1, ..., d
            if sum(1:d ~= uniqueLeaves) > 0
                error(['HACopula.checkleafs: All leaves of the HAC must be the set {1, ..., ' ...
                    num2str(d) '}.']);
            end
        end
        
        
        
        function checksnc(obj)
            %CHECKSNC - checks for the SNC in all parent-child pairs
            
            % check if obj.Model is a generator
            if ~isgenerator(obj.Family, obj.Parameter)
                error(['HACopula.checksnc: (' obj.Family ', ' num2str(obj.Parameter) ') is not a(n allowed) generator.']);
            end
            % check SNC with its children
            for iChild = 1:size(obj.Child,2)
                if isa(obj.Child{iChild},'HACopula')
                    % if the child is a fork, check the snc
                    checksncrec(obj, obj.Child{iChild});
                end
            end
        end
        
        
        
        function checksncrec(parent, child)
            %CHECKSNCREC - An auxiliary method for checksnc, which recursively checks the
            % SNC
            
            % is the child a generator
            if ~isgenerator(child.Family, child.Parameter)
                error(['HACopula.checksncrec: (' child.Family ', ' num2str(child.Parameter) ') is not a supported generator.']);
            end
            
            isSncSatisfied = true;
            
            % check for homogeneous vs heterogeneous case
            if strcmp(parent.Family, child.Family)
                % homogeneous case
                if strcmp(parent.Family, '14')
                    error(['HACopula.checksncrec: The sufficient nesting condition check for the family combination (' ...
                        parent.Family ', ' child.Family ') at levels ' num2str(parent.Level) ' and ' num2str(child.Level) ...
                        ' is not implemented (it is not known).']);
                else
                    if ~(parent.Parameter <= child.Parameter)
                        isSncSatisfied = false;
                    end
                end
                
            else
                % heterogeneous case
                switch [parent.Family child.Family]
                    case 'A19'                                          % the L_1 class
                        % is always satisfied
                    case {'C12', 'C19'}                                 % the L_2 class
                        if ~(parent.Parameter <= 1)
                            isSncSatisfied = false;
                        end
                    case {'AC', 'A20'}                                  % the L_3 class
                        if ~(child.Parameter >= 1)
                            isSncSatisfied = false;
                        end
                    case 'C14'
                        if ~(parent.Parameter * child.Parameter <= 1)   % the L_4 class
                            isSncSatisfied = false;
                        end
                    case 'C20'                                          % the L_4 class
                        if ~(parent.Parameter <= child.Parameter)
                            isSncSatisfied = false;
                        end
                    otherwise
                        error(['HACopula.checksncrec: The sufficient nesting condition check for the family combination (' ...
                            parent.Family ', ' child.Family ') at levels ' num2str(parent.Level) ' and ' num2str(child.Level) ...
                            ' is not implemented (the condition does not hold for any combination of the parameters or is not known).']);
                end
                
            end
            
            if ~isSncSatisfied
                error(['HACopula.checksncrec: The parent-child pair ((' parent.Family ', ' num2str(parent.Parameter) ...
                    '), (' child.Family ', ' num2str(child.Parameter) ...
                    ')) at levels ' num2str(parent.Level) ' and ' num2str(child.Level) ...
                    ' does not satisfy the sufficient nesting condition.']);
            end
            
            % check recursively SNC with the children of child
            for iChild = 1:size(child.Child,2)
                if isa(child.Child{iChild},'HACopula')
                    % if the child.Child is a fork, check the snc
                    checksncrec(child, child.Child{iChild});
                end
            end
        end
        
        
        
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
            % reEstimateType         - 'Ktauavg' - the [Górecki et al.,
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
            
            if ~(strcmp(reEstimateType, 'Ktauavg') || strcmp(reEstimateType, 'taumin'))
                error('HACopula::collapse: the parameter reEstimateType must be  either ''Ktauavg'' or ''taumin''.')
            end
            
            nForks = size(obj.Forks, 2);
            
            if nForks == 1
                error('HACopula::collapse: Attempting to collapse HAC that is AC.')
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
                    % use the re-estimation proposed in [Górecki et al., 2016b]
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
        
        
        
        function [areEqual, areEqualFamilies, nEqualFamilies] = comparestructures(HACObject1, HACObject2)
            %COMPARESTRUCTURES - returns areEqual = 1, if the structures of HACObject1 and
            % HACObject2 are equal, otherwise, returns areEqual = 0. 
            %
            % The comparison is based on comparing the sets of the
            % descendant leaves of all forks. Additionally, if the
            % structures are equal, the method compares the families of the
            % corresponding generators, and returns areEqualFamilies = 1,
            % if the families are all equal, areEqualFamilies = 0, at least
            % one family is different, nEqualFamilies is the number of the
            % equal families. If the structures are not equal,
            % areEqualFamilies = nEqualFamilies = NaN.
            
            d1 = getdimension(HACObject1);
            d2 = getdimension(HACObject2);

            areEqualFamilies = NaN;
            nEqualFamilies = NaN;            
            
            if (d1 ~= d2)
                areEqual = false;
                return;
            end
            
            d = d1;
            
            forks1 = HACObject1.Forks;
            forks2 = HACObject2.Forks;
            
            if (size(forks1,2) ~= size(forks2,2))
                areEqual = false;
                return;
            end
            
            k = size(forks1,2);
            
            % get leaves of each of the forks1
            leaves1 = cell(1, k);
            leaves1Extend = zeros(k, d);
            families1 = cell(1, k);
            leaves2 = cell(1, k);
            leaves2Extend = zeros(k, d);
            families2 = cell(1, k);
            for i = 1:size(forks1,2)
                leaves1{i} = forks1{i}.Leaves;
                leaves2{i} = forks2{i}.Leaves;
                
                % extend leaves for sortrow below
                leaves1Extend(i, 1:d) =  [zeros(1, d - size(leaves1{i},2)) sort(leaves1{i})];
                leaves2Extend(i, 1:d) =  [zeros(1, d - size(leaves2{i},2)) sort(leaves2{i})];
                
                % get the families
                families1{i} = forks1{i}.Family;
                families2{i} = forks2{i}.Family;
            end
            
            [leaves1Sorted, index1] = sortrows(leaves1Extend);
            [leaves2Sorted, index2] = sortrows(leaves2Extend);
            areEqual = (sum(sum(leaves1Sorted == leaves2Sorted)) == k*d);
            if areEqual
                sortedFamilies1 = families1(index1);
                sortedFamilies2 = families2(index2);
                nEqualFamilies = sum(strcmp(sortedFamilies1, sortedFamilies2));
                areEqualFamilies = (nEqualFamilies == k);
            else
               nEqualFamilies = NaN;
               areEqualFamilies = NaN;
            end
        end
        
        
        function UTrans = computediagonaltransform(obj, U)
            %COMPUTEDIAGONALTRANSFORM - transformes U(:, obj.Leaves) according to obj
            % into a single vector of observations through the diagonal
            % transformation (see [Górecki et al., 2014]).
            %
            % TODO: use this method in collapse for MLE in the diagonal
            % estimation 
            
            % prepare necessary functions
            psi = getgenerator(obj.Family, obj.Parameter, 0);
            psiInv = getgenerator(obj.Family, obj.Parameter, 1);
            nChild = length(obj.Child);
            delta = @(t) psi(nChild.* psiInv(t));  % the transformation
            
            % prepare necessary data
            UChild = zeros(size(U,1), nChild);
            for i = 1:nChild
                if isa(obj.Child{i}, 'HACopula')
                    % go into recursion
                    UChild(:, i) = computediagonaltransform(obj.Child{i}, U);
                else
                    % use the observations from U
                    UChild(:, i) = U(:, obj.Child{i});
                end
            end
            
            % compute the transformation
            UTrans = delta(max(UChild,[], 2)); 
            
            % do a NaN check
            [UTrans, nNaNs] = nanapprox(UTrans, UChild);
            if nNaNs > 0
                warning(['HACopula::computediagonaltransform: ' num2str(nNaNs) ' NaNs detected in the diagonal transformation and replaced by their approximations.']);
            end
        end
       
        
        
        function decreaselevel(obj)
            %DECREASELEVEL - decreases the level of all forks in obj. 
            %
            % Used after collapsing of two forks in the method collapse.
            %
            % NOTE: This method modifies the input obj.
            
            for i = 1:size(obj.Forks,2)
                if isa(obj.Forks{i},'HACopula')
                    obj.Forks{i}.Level = obj.Forks{i}.Level - 1;
                end
            end
        end
        
        
        
        function dist = distance(obj, objToCompare, varargin)
            %DISTANCE - computes a standardized Euclidian-like distance between a
            % selected type of matrix (kendall, upper-tail, lower-tail)
            % corresponding to obj and to objToCompare, where objToCompare
            % can be either an instance of the HACopula class or a Kendall
            % correlation matrix
            %
            % Note:
            % The type of matrix can be supplied as the third argument. If it
            % is not supplied, the default value is 'kendall'. Also, if
            % objToCompare is a Kendall correlation matrix, do not supply
            % the third argument as the 'kendall' type will be used.
            
            % Inputs checking
            narginchk(2,3);
            d = getdimension(obj);
            if isa(objToCompare,'HACopula')
                if d ~= getdimension(objToCompare)
                    error('HACopula::distance: Both HACopula inputs must be of the same dimension.');
                end
                type = varargin{1};
                if sum(strcmp({'kendall', 'upper-tail', 'lower-tail'}, type)) == 0
                    error('HACopula::distance: The input *type* must be from {''kendall'', ''upper-tail'',''lower-tail''}.');
                end 
            elseif ismatrix(objToCompare)
                if ~((max(size(objToCompare)) == min(size(objToCompare))) && ...
                        (max(size(objToCompare)) == d))
                    error('HACopula::distance: The matrix objToCompare must be a of size d*d.');
                end
                type = 'kendall';
                if ~isempty(varargin)
                    warning('HACopula::distance: Ignoring the third argument and setting type=''kendall''.')
                end
            else
                error('HACopula::distance: The input objToCompare must be an instance of HACopula class or a Kendall correlation matrix.');
            end
            
            % Distance computation
            dist = 0;
            for i = 1:d
                for j = i+1:d
                    [~, fork] = getyca(obj, i, j);
                    if isa(objToCompare,'HACopula')
                        [~, forkToCompare] = getyca(objToCompare, i, j);
                        switch type
                            case 'kendall'
                                dist = dist + (fork.Tau - forkToCompare.Tau)^2;
                            case 'upper-tail'
                                dist = dist + ...
                                       (gettaildependence(fork.Family, fork.Parameter, 'upper') - ...
                                        gettaildependence(forkToCompare.Family, forkToCompare.Parameter, 'upper'))^2;
                            case 'lower-tail'
                                dist = dist + ...
                                       (gettaildependence(fork.Family, fork.Parameter, 'lower') - ...
                                        gettaildependence(forkToCompare.Family, forkToCompare.Parameter, 'lower'))^2;
                        end
                    else % objToCompare is a matrix d*d
                        % type is kendall
                        dist = dist + (fork.Tau - objToCompare(i,j))^2;
                    end
                end
            end
            % standardization
            dist = sqrt(dist/nchoosek(d,2));
        end

        function out = evalsurv(obj, U)
        %EVALSURV - evaluates the survival copula of obj at U \in [0, 1]^d
        %
        % NOTE:
        % The computation involves evaluations of the HAC obj in all 2^d
        % corners of the d-dimensional hypercube (U, 1).
            
            out = prob(obj, U, 1*ones(1, getdimension(obj)));
        end
        
        function out = evaluate(obj, U)
            %EVALUATE - evaluates the HAC obj at a certain point (U)
            %
            % Assuming that obj represents a copula function C, the output 
            % is C(U), where C is applied row-wise on the rows of U.
            
            n = size(U, 1);
            nChildren = size(obj.Child,2);
            subU = zeros(n, nChildren);
            for i = 1:nChildren
                if isa(obj.Child{i}, 'HACopula')
                    % go to recursion
                    subU(:,i) = evaluate(obj.Child{i}, U);
                else % it is a leaf
                    % return the column corresponding to the leaf
                    subU(:,i) = U(:,obj.Child{i});
                end
            end
            psi = getgenerator(obj.Family, obj.Parameter, 0);
            psiInv = getgenerator(obj.Family, obj.Parameter, 1);
            out = psi(sum(psiInv(subU), 2));
        end
        
        
        
        function HACCdf = getcdf(obj)
            %HACCDF - a symbolic representation of the CDF
            
            syms t;
            syms theta;
            
            theta_i = sym(sprintf('theta%d',obj.TauOrdering));
            
            psi = getsymbgenerator(obj.Family, 0);      %get generator
            psi = subs(psi, theta, theta_i);
            psiInv = getsymbgenerator(obj.Family, 1);  %get generator inverse
            psiInv = subs(psiInv, theta, theta_i);
            
            psiInvSum = 0;
            for i = 1:length(obj.Child)
                
                if ~isa(obj.Child{i},'HACopula')
                    % the child is a leaf
                    child = sym(sprintf('u%d',obj.Child{i}));
                    psiInvSum = psiInvSum + subs(psiInv, t, child);
                else
                    child = getcdf(obj.Child{i});
                    psiInvSum = psiInvSum + subs(psiInv, t, child);
                end
            end
            
            HACCdf = subs(psi, t, psiInvSum);
        end
        
        
        
        function dimension = getdimension(obj)
            %GETDIMENSION - returns the dimension of the HAC obj.
            
            dimension = size(obj.Leaves,2);
        end
              
        
        
        function maxLevel = getmaxlevel(obj)
            %GETMAXLEVEL - the maximal nesting level in obj.
            
            maxLevel = -Inf;
            for i = 1:size(obj.Forks,2)
                maxLevel = max(maxLevel, obj.Forks{i}.Level);
            end
        end
        

        
        function families = getfamilies(obj)
            %GETFAMILIES - the set of families involved in the obj.
            
            families = {};
            for i = 1:length(obj.Forks)
                family = obj.Forks{i}.Family;
                iFamily = find(strncmp(families, family, 2)); % is the family already in families
                if size(iFamily,2) == 0
                    % add this family
                    families = [{family} families];
                end
            end
        end
        
        
        
        function fork = getforkwithtauindex(obj, tauIndex)
            %GETFORKWITHTAUINDEX - returns the fork with TauOrdering equal to tauIndex.
            
            for i = 1:length(obj.Forks)
                if obj.Forks{i}.TauOrdering == tauIndex
                    fork = obj.Forks{i};
                    return
                end
            end
        end
        
        
        
        function HACPdf = getpdf(obj)
            %HACPFD - a symbolic representation of the PDF.
            
            leaves = obj.Leaves;
            HACCdf = getcdf(obj);
            %differentiate HACCdf
            HACPdf = HACCdf;
            for i = 1:length(leaves)
                HACPdf = diff(HACPdf, sym(sprintf('u%d',leaves(i))));
            end
        end
        
        
        
        function depMatrix = getdependencematrix(obj, type)
            %GETDEPENDENCEMATRIX - returns a [d x d] matrix of a selected dependence
            % coefficients (according to type from {kendall, lower-tail and
            % upper-tail}). If type=tails, then tail coefficients are
            % condensed into one matrix and the upper-tail coefficients are
            % placed above the diagonal, whereas the lower-tail
            % coefficients below it.
            %
            % NOTE:
            % The main diagonal always contains ones as all the considered
            % coefficients are 1 for two same random variables.
            
            if sum(strcmp({'kendall', 'upper-tail', 'lower-tail', 'tails'}, type)) == 0
                error('HACopula::getdependencematrix: The input *type* must be from {''kendall'', ''upper-tail'', ''lower-tail'', ''tails''}.')
            end
            d = getdimension(obj);
            depMatrix = ones(d, d);
            for i = 1:d
                for j = i+1:d
                    [~, fork] = getyca(obj, i, j);
                    switch type
                        case 'kendall'
                            depMatrix(i, j) = fork.Tau;
                            depMatrix(j, i) = depMatrix(i, j);
                        case 'upper-tail'
                            depMatrix(i, j) = gettaildependence(fork.Family, fork.Parameter, 'upper');
                            depMatrix(j, i) = depMatrix(i, j);
                        case 'lower-tail'
                            depMatrix(i, j) = gettaildependence(fork.Family, fork.Parameter, 'lower');
                            depMatrix(j, i) = depMatrix(i, j);
                        case 'tails'
                            depMatrix(i, j) = gettaildependence(fork.Family, fork.Parameter, 'upper');
                            depMatrix(j, i) = gettaildependence(fork.Family, fork.Parameter, 'lower');
                    end
                end
            end
        end
        
       
        function [ycaTauOrdering, ycaFork] = getyca(obj, leaf1, leaf2)
            %GETYCA - returns as ycaTauOrdering the TauOrdering property of the youngest common
            % ancestor of the leaves leaf1 and leaf2.
            % Also the youngest common ancestor fork is returned in ycaFork
            
            if (size(intersect([leaf1 leaf2], obj.Leaves), 2) ~= 2)
                error('HACopula::getycaindex: leaf1 and leaf2 must be two different leaves of the HAC.');
            end
            orderedHACModel = obj.Forks;
            d = getdimension(obj);
            k = size(orderedHACModel, 2);
            minLeavesSize = d+1;
            % get all sets of the leaves
            for i = 1:k
                leavesVec = orderedHACModel{i}.Leaves;
                if (size(intersect([leaf1 leaf2], leavesVec), 2) == 2) && ...
                        (size(leavesVec, 2) < minLeavesSize)
                    % yca is the fork with the fewest descendat leaves
                    % containing both leaf1 and leaf2
                    ycaTauOrdering = orderedHACModel{i}.TauOrdering;
                    ycaFork = orderedHACModel{i};
                end
            end
        end
        
        
        
        function Sn = goffullpairwise(obj, U, aggFcn, statisticName, emp2copulas)
            %GOFFULLPAIRWISE - Aggregated pairwise goodness-of-fit
            %statistic
            %
            % This function computes the bivariate version of the Cramér-von Mises 
            % statistic statisticName ('K', 'R' or 'E') for each of the
            % bivariate margins of the HAC obj and the observations U, and
            % aggregates them using the aggregation function aggFcn
            % ('average', 'max' or 'mean'). 
            % 
            % NOTE:
            % The input emp2copulas is optional and can be pre-computed using
            % computeallemp2copulas(U).
            
            % perform basic data checks
            if ~iscopuladata(U)
                warning('HACopulafit::goffullpairwise: At least one column of U was rejected to be standard uniform according to KS test at the 5% significance level, i.e., U might not be a sample from a copula.');
            end
            
            if (~exist('emp2copulas', 'var'))
                emp2copulas = computeallemp2copulas(U);
            end
            
            d = size(U, 2);
            if getdimension(obj) ~= d
                error('HACopula::goffullpairwise: The data and the HAC are of different dimensions.')
            end
            
            pairSn = zeros(d, d);
            for i = 1:d
                for j = i+1:d
                    fork = getforkwithtauindex(obj, getyca(obj, i, j));
                    pairSn(i, j) = gof2Sn(U(:, [i j]), fork.Family, fork.Parameter, statisticName, emp2copulas{i,j});
                end
            end
            
            % get upper triangular part of pairSn (without the main diagonal)
            B = pairSn';
            pairSnVec = B(tril(true(size(B)),-1));
            
            switch aggFcn
                case 'max'
                    Sn = max(pairSnVec);
                case 'average'
                    Sn = mean(pairSnVec);
                case 'min'
                    Sn = min(pairSnVec);
                otherwise
                    error('gof2Sn: Unsupported aggregation function.');
            end
        end
        
        
        
        function Sn = gofdSnE(obj, U)
            %GOFDSNE - returns the value of the Cramér-von Mises statistics (denoted
            % S_n in [Genest et al., 2009]) based on the empirical copula for
            % the HAC obj and the observations U.
            
            % perform basic data checks
            if ~iscopuladata(U)
                warning('HACopulafit::gofdSnE: At least one column of U was rejected to be standard uniform according to KS test at the 5% significance level, i.e., U might not be a sample from a copula.');
            end
            
            %compute empirical copula values in data point
            [n, d] = size(U);
            
            Cn = zeros(n,1);
            for i=1:n
                mult = (U(:,1) <= U(i,1));
                for j=2:d
                    mult = mult .* (U(:,j) <= U(i,j));
                end
                Cn(i) = sum(mult);
            end
            Cn = Cn / n;
            
            % slow symbolic evaluation
            %get theoretical cdf of the copula for estimated parameter
            %             cdfCtheta = getcdf(obj);
            
            %compute statistics
            %             forksArray = getforksarray(obj);
            %             Usymb = sym('u%d', [1 d]);
            %             theta = sym('theta%d', [1 2*d-1]);
            %Ctheta = evalsymbformula(cdfCtheta, [Usymb, theta([forksArray(:).TauOrdering])], U, [forksArray(:).Parameter]);
            
            % fast non-symbolic evaluation
            Ctheta = evaluate(obj, U);
            
            % NaN check
            [Ctheta, nNaNs] = nanapprox(Ctheta, U);
            if nNaNs > 0
                warning(['HACopula:gofdSnE:: ' num2str(nNaNs) ' NaNs detected and replaced by their approximations.']);
            end
            
            % Inf check
            infs = isinf(Ctheta);
            nInfs = sum(infs);
            if (nInfs > 0)
                warning(['HACopula:gofdSnE:: ' num2str(nInfs) ' Infs detected and each is replaced by 1.']);
                % replace Inf by 1
                Ctheta(infs) = 1;
            end
            
            Sn = sum((Ctheta - Cn).^2);
        end
        
        
        
        function fig = plot(obj, varargin)
            %PLOT - Visualize a HAC.
            % 
            % Purpose:
            % Plots a tree- or a dendrogram-based representation of
            % a HACopula object obj.
            %
            % NOTE:
            % 1) The children nodes of each fork are orderer (left to right)
            % using the function setlrordering.
            % 2) The figure is positioned in the middle of the screen and
            % sized 1/2 of the width and height of the screen.
            % 3) The figure is plotted in two phases, where the edges are
            % plotted during the first phase and the nodes during the second.
            % 4) If needed, the description of the forks can be customized
            % through genText, whereas the description of the leaves can be
            % customized through varText in the auxiliary function
            % plotHACrec. Also, the size of the font and of the nodes can be
            % customized there.
            %
            % Inputs (obligatory):
            % obj               - a HACopula object
            % Inputs (optional):
            % GraphType         - A type of the graph ('tree' (default) or
            %                     'dendrogram').
            % Positioning       - A type of vertical positioning of the forks. 
            %                     If set to 'even' (default), the forks are
            %                     vertically positioned evenly according
            %                     to their Level property.
            %                     If set to 'tau, the forks are
            %                     vertically positioned according
            %                     to their Tau property.
            % ForkBackground    - If set to 'white' (default), each fork has a white
            %                     bar background. If set to 'none', no
            %                     background is drawn for a fork.
            % MarkerSize        - a positive integer determining the size of the
            %                     circle markers used for the leaves
            % FontSize          - a positive integer determining the font
            %                     size
            
            % default setting of the optional inputs 
            graphType = 'tree';
            positioning = 'even';
            forkBackground = 'white';
            markerSize = 35; % use 32 for d = 20, 35 otherwise
            fontSize = 17; % use 12 for d = 20, 20 otherwise

            narginchk(1,7);

            % check for additional parameters
            if mod(size(varargin,2),2) == 1
                error('HACopula.plot: there must be an even number of the additional parameters');
            end
            
            if size(varargin,2) >= 2
                % there are some additional parameters
                parNames = lower(varargin(1:2:size(varargin,2)-1));
                parValues = varargin(2:2:size(varargin,2));
                
                % check of the parameter names are allowed
                for i = 1:size(parNames,2)
                    try
                        validatestring(parNames{i},{'GraphType', 'Positioning', 'ForkBackground', 'MarkerSize', 'FontSize'});
                    catch
                        error(['The input, ''' parNames{i} ''', did not match any of the valid parameter names (GraphType, Positioning, ForkBackground, MarkerSize, FontSize).']);
                    end
                end
                
                % check for GraphType parameter
                iGraphType = find(strncmp(parNames, lower('GraphType'), 9));
                
                % is GraphType a parameter ?
                if size(iGraphType,2) > 0 %
                    if size(iGraphType,2) > 1
                        error('HACopula.plot: GraphType is a repeating parameter.')
                    end
                    graphType = parValues{iGraphType};
                    if ~(strcmp(graphType, 'tree') || strcmp(graphType, 'dendrogram'))
                        error('HACopula.plot: the value corresponding to GraphType must be ''tree'' or ''dendrogram''.');
                    end
                end
                
                % check for Positioning parameter
                iPositioning = find(strncmp(parNames, lower('Positioning'), 11));
                
                % is Positioning a parameter ?
                if size(iPositioning,2) > 0 %
                    if size(iPositioning,2) > 1
                        error('HACopula.plot: Positioning is a repeating parameter.')
                    end
                    positioning = parValues{iPositioning};
                    if ~(strcmp(positioning, 'even') || strcmp(positioning, 'tau'))
                        error('HACopula.plot: the value corresponding to Positioning must be ''even'' or ''tau''.');
                    end
                    if strcmp(graphType, 'tree') && strcmp(positioning, 'tau')
                        warning(['HACopula.plot: Some of the edges could cross each other for the setting GraphType=''tree'' '...
                            'and Positioning=''tau''. Use GraphType=''tree'' and Positioning=''even'' instead or use GraphType=''dendrogram''.']);
                    end
                end
                
                % check for ForkBackground parameter
                iForkBackground = find(strncmp(parNames, lower('ForkBackground'), 13));
                
                % is ForkBackground a parameter ?
                if size(iForkBackground,2) > 0 %
                    if size(iForkBackground,2) > 1
                        error('HACopula.plot: ForkBackground is a repeating parameter.')
                    end
                    forkBackground = parValues{iForkBackground};
                    if ~(strcmp(forkBackground, 'white') || strcmp(forkBackground, 'none'))
                        error('HACopula.plot: the value corresponding to ForkBackground must be ''white'' or ''none''.');
                    end
                end
                
                % check for MarkerSize parameter
                iMarkerSize = find(strncmp(parNames, lower('MarkerSize'), 10));                
                % is MarkerSize a parameter ?
                if size(iMarkerSize,2) > 0 %
                    if size(iMarkerSize,2) > 1
                        error('HACopula.plot: MarkerSize is a repeating parameter.')
                    end
                    markerSize = parValues{iMarkerSize};
                end
 
                % check for FontSize parameter
                iFontSize = find(strncmp(parNames, lower('FontSize'), 10));                
                % is FontSize a parameter ?
                if size(iFontSize,2) > 0 %
                    if size(iFontSize,2) > 1
                        error('HACopula.plot: FontSize is a repeating parameter.')
                    end
                    fontSize = parValues{iFontSize};
                end
                
                
            end
            
            
            % get models parameters
            maxLevel = getmaxlevel(obj);
            d = getdimension(obj);
            k = size(obj.Forks,2);
            
            % do a copy of HACopula in order not to change the original
            % object by the function setlrordered
            cpObj = copy(obj);
            % interchange children according to numbers of descendant forks
            setlrordering(cpObj);
            
            % start plotting
            fig = figure;
            cla reset;
            screensize = get( 0, 'Screensize' );
            set(fig, 'Position', [round(screensize(3)/4), round(screensize(4)/4),...
                round(screensize(3)/2), round(screensize(4)/2)]);
            set(gcf,'Units','normal')
            %set(gca,'Position',[0.1 0.08 0.88 0.85])
            xMarg = 95/screensize(3);
            yMarg = 100/screensize(4);
            set(gca, 'Position', [xMarg/2 yMarg 1-2*xMarg 1-2*yMarg]);
            %set(gca,'ButtonDownFcn','selectmoveresize');
            
            hold on;
            
            % do the recursive plot
            plotHACrec(cpObj, 0, maxLevel, 1, d, k, positioning, graphType, forkBackground, 0, markerSize, fontSize);
            plotHACrec(cpObj, 0, maxLevel, 0, d, k, positioning, graphType, forkBackground, 0, markerSize, fontSize);
            axis off;
            
            %format figure
            switch positioning
                case 'even'
                    % add a vertical axis showing the levels in a HAC
%                     axis([0.98 d (-maxLevel - 1.1) -0.7 ]); %[xmin xmax ymin ymax]
%                     move = 0.1; % shifts horizontally the vertical axis
%                     set(gca,'XLim',[move - 0.05 d]);
%                     
%                     plot([move move], [-1 -maxLevel], 'k');
%                     for i = 1:maxLevel
%                         plot([move (-0.05 + move)], [-i -i], 'k');
%                         text(move-0.05, -i, num2str(i), 'HorizontalAlignment', 'right', 'FontSize',17, 'Interpreter','latex', 'rotation',0);
%                     end
%                     h = text(move-0.5, -(maxLevel-1)/2-1, 'Level', 'HorizontalAlignment', 'center', 'FontSize',17, 'Interpreter','latex', 'rotation',0);
%                     set(h, 'rotation', 90)
                case 'tau'
                    axis([0.98 d 0 1]); %[xmin xmax ymin ymax]
                    move = 0.3;
                    %set(gca,'YLim',[0 1],'Layer','top');
                    set(gca,'XLim',[move - 0.05 d]);
                    
                    plot([move move], [0 1], 'k');
                    plot([move (-0.05 + move)], [0 0], 'k');
                    plot([move (-0.05 + move)], [1 1], 'k');
                    text(move, 0.5, '$\tau$ ', 'HorizontalAlignment', 'right', 'FontSize',17, 'Interpreter','latex', 'rotation',0);
                    text(move, 0, '1 ', 'HorizontalAlignment', 'right', 'FontSize',17, 'Interpreter','latex', 'rotation',0);
                    text(move, 1, '0 ', 'HorizontalAlignment', 'right', 'FontSize',17, 'Interpreter','latex', 'rotation',0);
                otherwise
                    error(['HACopula::plotHACrec: Positioning ' positioning ' is not supported.']);
            end
            
            hold off;
            
            % save the plot to a specified location
            %set(fig_handle1, 'PaperPositionMode','auto');
            %print(fig_handle1, '-depsc', ['k:\copulas\figures\copula_' ...
            %                               num2str(d) 'D_' [ord(:).type] '.eps']);
        end
        
        function out = prob(obj, l, u)
            % PROP - probability of falling into a hypercube from [0,1]^d
            %
            % Compute the probability of a d-dimensional random vector U
            %   distributed according to a given copula *obj* to fall in a
            %   hypercube (l, u], where *l* and *u* denote the lower and upper
            %   corners of the hypercube, respectively.
            %
            % NOTE:
            % The computation involves evaluations of the HAC obj in all 2^d
            % corners of the hypercube.
            
            d = getdimension(obj);
            % get the corners of the d-dimensional unit hypercube.
            [corner{1:d}] = ndgrid(logical([0 1]));
            corner = cat(d+1,corner{d:-1:1});
            corner = double(reshape(corner,[],d));
            
            % get the corners of [l u]^d
            lu = zeros(size(corner));
            for i = 1:d
                lu(:, i) = corner(:, i) * (u(i)-l(i)) + l(i);
            end
            
            % evaluate the copula at these corners
            cornProb = evaluate(obj, lu);
            
            % compute the C-volume of [l u]^d
            out = 0;
            for i = 1:size(corner, 1)
                out = out + (-1)^sum(1-corner(i,:)) * cornProb(i);
            end
        end
        
        
        function observations = rnd(obj, n)
            %RND - Sample n observations (rows) according to obj.
            
            observations = HACopularnd(obj, n);
        end
        
        
        
        function pvalue = computepvalue(obj, U, estimator, N)
            % COMPUTEPVALUE - bootstrap computation of p-value
            %
            % Purpose:
            % Given a HAC estimate, pseoudo-observations and an estimator,
            % it performs the bootstrap method proposed in [Genest et al.,
            % 2009] in order to get a p-value estimate.
            %
            % Inputs:
            % obj       - A HAC estimate obtained by the input estimator
            %             for the input data U.
            % U         - Pseudo-observations.
            % estimator - A HAC estimator with only one input argument 
            %             corresponding to pseudoobservations, e.g., an
            %             anonymous funtion @(U)f(U).
            % N         - A desired number of bootstrap replications.
            %
            % Outputs:
            % pvalue    - A value from [0, 1], which can be interpreted as 
            %             the probability (under valid H0, i.e., the
            %             underlying copula is obj) that, using a given HAC
            %             estimator, the goodness-of-fit, would be worse
            %             than the actual observed results, i.e.,
            %             goodness-of-fit of the
            %             estimate obj (obtained by the input estimator)
            %             for the data U.
            %
            % NOTE:
            % The goodness-of-fit statistic used in the computation is the
            % statistic denoted S_n in [Genest et al., 2009].
            
            n = size(U,1);
            estSn = gofdSnE(obj, U);
            [~,theirSn] = bootstrapSn(obj, estimator, N, n);
            pvalue = sum(theirSn > estSn)/size(theirSn,2);
        end
        
        
        function setlrordering(obj)
            %SETLRORDERING - switch children of a fork in order to canonize
            %                the graphical representation of a HAC. 
            %
            % This auxiliary function (for the plot function) interchanges
            % the children of each fork in order to get a left-right
            % ordering, i.e., a child with less dencendant forks are placed
            % at the left-side of a child with more dencendant forks. If
            % the numbers of the dencendant forks are equal for two or more
            % children, the child with the lowest descendant leaf index is
            % placed at the left-side of the remaining children (which are
            % ordered analogously).
            %
            % This ordering assures that a fully-nested Archimedean copula
            % is plotted (by the plot function) as is common on "the main
            % diagonal"
            %    
            % NOTE: This function modifies the input obj.
            
            nChildren = size(obj.Child, 2);
            nDescForks = zeros(1, nChildren);
            lowestLeafIndex = zeros(1, nChildren);
            for i = 1:nChildren
                if isa(obj.Child{i},'HACopula')
                    nDescForks(i) = size(obj.Child{i}.Forks,2) - 1;
                    lowestLeafIndex(i) = min(obj.Child{i}.Leaves);
                else
                    nDescForks(i) = -1;
                    lowestLeafIndex(i) = obj.Child{i};
                end
            end
            [~, lrOrdering] = sortrows([nDescForks' lowestLeafIndex']);
            % interchange the children according to lrOrdering
            child = cell(1, nChildren);
            for i = 1:nChildren
                child{i} = obj.Child{lrOrdering(i)};
            end
            obj.Child = child;
            
            % do recursion
            for i = 1:nChildren
                if isa(obj.Child{i},'HACopula')
                    setlrordering(obj.Child{i});
                end
            end
        end
        
        
        
        function out = tolatex(obj, type)
            %TOLATEX - returns a latex expression of cdf or pdf of the HAC in *obj*,
            % i.e., the input *type* should be in {'cdf', 'pdf'}.
            %
            % NOTE:
            % Be aware that for d > 5, getting pdf becomes extremely time
            % demanding
            
            if strcmp(type, 'cdf')
                f = getcdf(obj);
            elseif strcmp(type, 'pdf')
                f = getpdf(obj);
            else
                error('HACopula: Unsupported type. Choose one from {''cdf'', ''pdf''}.');
            end
            
            %out = latex(simplify(f)); % simplification - uncommnent to simplify the formula
            %out = latex(simplify(f,'IgnoreAnalyticConstraints',true)); % simplification - uncommnent to simplify the formula
            out = latex(f);

            d = getdimension(obj);
            for i = 1:d
                out = regexprep(out, sprintf('%s%d%s','\\mathrm{u', i, '}'), sprintf('%s%d%s','u_{', i, '}'));
            end
            for i = d+1:d+length(obj.Forks)
                out = regexprep(out, sprintf('%s%d%s','\\mathrm{theta', i, '}'), sprintf('%s%d%s','\\theta_{', i, '}'));
            end
        end

        
        
    end
end



% -------------------------------------------------------------------------
% end of the HACopula classdef




% -------------------------------------------------------------------------
% auxiliary functions
% -------------------------------------------------------------------------


function [out1, out2] = plotHACrec(HACModel, iVariable, maxLevel, lines, d, k,...
    positioning, graphType, forkBackground, parentLevel, markerSize, fontSize)
% An auxiliary recursive function for the HACopula's method plot.


% customize the look of the plot
STEP_X = 1; % a horizontal distance of two nodes next to each other
STEP_Y = -1; % % a vertical distance of two nodes from levels i and i+1, respectively

if isa(HACModel, 'HACopula')
    nChildren = size(HACModel.Child,2);
else
    nChildren = 0;
end

%plot the edges
if (nChildren > 1)
    x = zeros(1, nChildren);
    for i = 1:nChildren
        [iVariable, x(i)] = plotHACrec(HACModel.Child{i}, iVariable, maxLevel, lines, d, k, positioning, graphType, forkBackground, HACModel.Level, markerSize, fontSize);
    end
    xm = mean([x(1) x(nChildren)]);
    
    if (lines == 1)
        for i = 1:nChildren
            if isnumeric(HACModel.Child{i})
                % it's a leaf
                switch positioning
                    case 'even'
                        switch graphType
                            case 'dendrogram'
                                plot([x(i) x(i) xm], [(maxLevel + 1) * STEP_Y  HACModel.Level ...
                                    * STEP_Y HACModel.Level * STEP_Y], 'k');
                            case 'tree'
                                plot([x(i) xm], [(HACModel.Level + 1) * STEP_Y  ...
                                    HACModel.Level * STEP_Y], 'k');
                            otherwise
                                error(['HACopula.plotHACrec: GraphType ' graphType ' is not supported.']);
                        end
                    case 'tau'
                        switch graphType
                            case 'dendrogram'
                                plot([x(i) x(i) xm], [0  1-HACModel.Tau ...
                                    1-HACModel.Tau], 'k');
                            case 'tree'
                                plot([x(i) xm], [0  1-HACModel.Tau], 'k');
                            otherwise
                                error(['HACopula.plotHACrec: GraphType ' graphType ' is not supported.']);
                        end
                    otherwise
                        error(['HACopula.plotHACrec: positioning ' positioning ' is not supported.']);
                end
            else
                % it is a fork
                switch positioning
                    case 'even'
                        switch graphType
                            case 'dendrogram'
                                plot([x(i) x(i) xm], [HACModel.Child{i}.Level * STEP_Y ...
                                    HACModel.Level * STEP_Y HACModel.Level * STEP_Y], 'k');
                            case 'tree'
                                plot([x(i) xm], [HACModel.Child{i}.Level * STEP_Y ...
                                    HACModel.Level * STEP_Y], 'k');
                            otherwise
                                error(['HACopula.plotHACrec: GraphType ' graphType ' is not supported.']);
                        end
                    case 'tau'
                        switch graphType
                            case 'dendrogram'
                                plot([x(i) x(i) xm], [1-HACModel.Child{i}.Tau  1-HACModel.Tau ...
                                    1-HACModel.Tau], 'k');
                            case 'tree'
                                plot([x(i) xm], [1-HACModel.Child{i}.Tau 1-HACModel.Tau], 'k');
                            otherwise
                                error(['HACopula.plotHACrec: GraphType ' graphType ' is not supported.']);
                        end
                    otherwise
                        error(['HACopula.plotHACrec: positioning ' positioning ' is not supported.']);
                end
            end
        end
    end
end

%plot the nodes
if isnumeric(HACModel) % is a leaf
    iVariable = iVariable + 1; %variable counter
    x = iVariable * STEP_X;
    switch positioning
        case 'even'
            switch graphType
                case 'dendrogram'
                    y = (maxLevel + 1) * STEP_Y;
                case 'tree'
                    y = (parentLevel + 1) * STEP_Y;
                otherwise
                    error(['HACopula.plotHACrec: GraphType ' graphType ' is not supported.']);
            end
        case 'tau'
            y = 0;
        otherwise
            error(['HACopula.plotHACrec: positioning ' positioning ' is not supported.']);
    end
    
    if (lines == 0)
        pl = plot(x, y, 'o');
        set(pl,'MarkerSize',markerSize, 'MarkerEdgeColor','k', 'MarkerFaceColor','w');
    end
    %plot_circle(x, y, marker_shift);
else
    x = xm;
    switch positioning
        case 'even'
            y = HACModel.Level * STEP_Y;
        case 'tau'
            y = 1-HACModel.Tau;
        otherwise
            error(['HACopula.plotHACrec: positioning ' positioning ' is not supported.']);
    end
end

if (lines == 0)
    %plot text
    if isnumeric(HACModel)
        type = 'u';
        add = '';
    else
        type = HACModel.Family;
        if type == '?'
            add = '?';
            whiteSpUnknown = '~~~';
        else
            add =  sprintf('%3.3f',HACModel.Parameter);
            whiteSpUnknown = '';
        end
        %add =  ['$\theta_{' num2str(HACModel.TauOrdering) '}$'];
    end
    
    if isnumeric(HACModel)
        varText = ['$' type '_{', num2str(HACModel), '}$'];
        text(x, y, varText, 'HorizontalAlignment', 'center', 'FontSize',fontSize, 'Interpreter','latex', 'FontUnits','pixels');
    else
        
        % do centering :)
        if HACModel.TauOrdering < 10
            whiteSpace = '~~~~';
        else
            whiteSpace = '~~~';
        end
        % generators (forks) representation
        % for the outputs for [Górecki et al., 2016b]
        if strcmp(HACModel.Family,'A') || strcmp(HACModel.Family,'C')
            whiteSpaceAdd = '\,';
        else
            whiteSpaceAdd = '~';
        end
        genText = ['$' whiteSpace '\lambda(' num2str(HACModel.TauOrdering) ')$' char(10) whiteSpUnknown '(' type ', ' add ')' char(10) whiteSpaceAdd '$\tau = ' sprintf('%3.3f',HACModel.Tau) '$'];  % [Górecki et al., 2016b] example

        % plot text
        text(x, y, genText,...
            'HorizontalAlignment', 'center', 'FontSize',fontSize, 'Interpreter','latex', 'FontUnits','pixels',...
            'BackgroundColor', forkBackground, 'Margin', 1);
    end
end

out1 = iVariable;
out2 = x;
end




