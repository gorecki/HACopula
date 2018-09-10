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
%
%
% Copyright 2018 Jan Gorecki

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
                error('HACopula:cell2tree', ['HACopula.cell2tree: The generator (' HACCellModel{i}{1} ', ' num2str(HACCellModel{i}{2}) ') is not supported. For the supported generators, see the function isgenerator.']);
            end
        end
    end
    if ~isFound
        error('HACopula:cell2tree', 'HACopula.cell2tree: wrong cell HAC structure.');
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
        
        obj.Dim = size(obj.Leaves,2);
    end
else
    error('HACopula:cell2tree', 'HACopula.cell2tree: Wrong cell HAC structure. Invalid number of subcells.')
end
end
