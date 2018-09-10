%classdef HACopula < matlab.mixin.Copyable
classdef HACopula < handle
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
    % e.g., see [Gorecki et al., 2017].
    %
    %
    % References:
    % [Genest et al., 2009] Genest, C., Remillard, B., and Beaudoin, D.
    %    (2009). Goodness-of-fit tests for copulas: A review and a power
    %    study. Insurance: Mathematics and Economics, 44(2):199-213.
    % [Gorecki et al., 2014] Gorecki, J., Hofert, M., and Holena, M. (2014). On
    %     the consistency of an estimator for hierarchical Archimedean copulas.
    %     In 32nd International Conference on Mathematical Methods in
    %     Economics, pages 239-244.
    % [Gorecki et al., 2017] On Structure, Family and Parameter Estimation
    %     of Hierarchical Archimedean copulas. Journal of Statistical Computation 
    %     and Simulation, 87(17), 3261-3324
    % [Nelsen, 2006] Nelsen, R. (2006). An Introduction to Copulas. Springer,
    %     2nd edition.
    % [Uyttendaele, et al., 2016] Uyttendaele, et al., On the estimation of
    %     nested Archimedean copulas: A theoretical and an experimental
    %     comparison, Tech. rep., UCL (2016).
    %
    %
    % Copyright 2018 Jan Gorecki
    
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
        Dim             % The dimension of the copula (which equals to the number of its leaves)
        % HACopula handles
        Child           % A cell array of child (HACopula) objects.
        Parent          % [], if Level = 1, otherwise it contains the (HACopula) parent of this object.
        Root            % The (HACopula) root of this  object.
        Forks           % A cell array with HACopula objects containing all 
                        % descendant forks of this fork including this
                        % fork.
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
                error('HACopula:constructor', 'HACopula: Unrecognized HAC model.');
            end
            % do basic checks
            checkleaves(obj);
            checksnc(obj);
            % add tau ordering
            addtauordering(obj);
        end
        
    end
end

