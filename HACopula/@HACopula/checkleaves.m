function checkleaves(obj)
%CHECKLEAVES - check if the leaves of a HACopula object consitute {1, ..., d}.
%
%
% Copyright 2018 Jan Gorecki

%get leaves of the model
leaves = obj.Leaves;

% check if the leaves are integers
if sum(floor(leaves) == leaves) ~= size(leaves,2)
    error('HACopula:checkleafs', 'HACopula.checkleafs: The leaves of the HAC are not integer numbers.');
end

uniqueLeaves = unique(leaves);
if size(uniqueLeaves,2) ~= size(leaves,2)
    error('HACopula:checkleafs', 'HACopula.checkleafs: The leaves of the HAC are not unique, i.e., there are two or more of the same leaves.');
end

% the dimension of the HAC is the number of its leaves
d = size(leaves, 2);

% the leaves must be 1, ..., d
if sum(1:d ~= uniqueLeaves) > 0
    error('HACopula:checkleafs', ['HACopula.checkleafs: All leaves of the HAC must be the set {1, ..., ' ...
        num2str(d) '}.']);
end
end
