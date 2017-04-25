function checksncrec(parent, child)
%CHECKSNCREC - An auxiliary method for checksnc, which recursively checks the
% SNC
%
%
% Copyright 2017 Jan Górecki

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