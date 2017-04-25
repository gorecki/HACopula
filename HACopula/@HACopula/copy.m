function cpObj = copy(obj)
% returns a copy of a HACopula object
%
%
% Copyright 2017 Jan Górecki

% Make a shallow copy of all properties
cpObj = HACopula();
cpObj.Family = obj.Family;
cpObj.Parameter = obj.Parameter;
cpObj.Tau = obj.Tau;
cpObj.TauOrdering = obj.TauOrdering;
cpObj.Level = obj.Level;
cpObj.Leaves = obj.Leaves;

% redirect to the copied forks
cpObj.Forks = {cpObj};
% Make a deep copy
for i = 1:size(obj.Child,2)
    if isa(obj.Child{i}, 'HACopula')
        cpObj.Child{i} = copy(obj.Child{i});
        % redirect to the copied parent
        cpObj.Child{i}.Parent = cpObj;
        % redirect to the copied forks
        cpObj.Forks = [cpObj.Forks cpObj.Child{i}.Forks];
    else
        cpObj.Child{i} = obj.Child{i}; % for leaves, no deep copy is needed
    end
end
% redirect to the copied root
if cpObj.Level == 1
    for i = 1:size(cpObj.Forks, 2)
        cpObj.Forks{i}.Root = cpObj;
    end
end
end