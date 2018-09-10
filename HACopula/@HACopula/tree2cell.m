function cellModel = tree2cell(obj)
%TREE2CELL - convert a HACopula object to a nested cell structure 
%
% This function performs a process that is inversion of what happens in
% cell2tree.
%
% Example:
% myCell = {{'A', 0.5}, 1, {{'C', 1.25}, 2, {{'20', 1.5}, 3, 4}, {{'C', 1.5}, 5, 6}}};
% myHAC = HACopula(myCell);
% myCell2 = tree2cell(myHAC);  
%
% now myCell2 and myCell are equal 


cellModel = cell(1, 1 + length(obj.Child));
cellModel{1} = {obj.Family, obj.Parameter};

for i = 1:length(obj.Child)
    if isa(obj.Child{i}, 'HACopula')
        cellModel{i + 1} = tree2cell(obj.Child{i});
    else % variable
        cellModel{i + 1} = obj.Child{i};
    end
end

end