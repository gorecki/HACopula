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
%
%
% Copyright 2018 Jan Gorecki

d = obj.Dim;

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
