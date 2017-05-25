function disp(obj)
% In Octave, emulates the output of MATLAB's function disp.
% In MATLAB, just uses the disp function provided by the superclass *handle*
%
%
% Copyright 2017 Jan Górecki

if isoctave
    disp(['  HACopula with properties:' char(10)]);
    
    childStr = '{';
    for i = 1:length(obj.Child)
        if isa(obj.Child{i}, 'HACopula')
            childStr = [childStr '[1x1 HACopula]'];
        else
            childStr = [childStr '[' num2str(obj.Child{i}) ']'];
        end
        if i < length(obj.Child)
            childStr = [childStr ' '];
        end
    end
    childStr = [childStr '}'];
    
    if isa(obj.Parent, 'HACopula')
        parentStr = '[1x1 HACopula]';
    else
        parentStr = '[]';
    end
    
    forksStr = '{';
    for i = 1:length(obj.Forks)
        forksStr = [forksStr '[1x1 HACopula]'];
        if i < length(obj.Forks)
            forksStr = [forksStr ' '];
        end
    end
    forksStr = [forksStr '}'];
    
    disp(['      Family: ' obj.Family char(10) ...
          '   Parameter: ' sprintf('%.4f', obj.Parameter) char(10) ...
          '         Tau: ' sprintf('%.4f', obj.Tau) char(10) ...
          ' TauOrdering: ' num2str(obj.TauOrdering) char(10) ...
          '       Level: ' num2str(obj.Level) char(10) ...
          '      Leaves: [' strjoin(strsplit(num2str(num2str(obj.Leaves)), '  '), ' ')  ']' char(10) ...
          '       Child: ' childStr char(10) ...
          '      Parent: ' parentStr char(10) ...
          '        Root: [1x1 HACopula]' char(10) ...
          '       Forks: ' forksStr char(10) ...
    ]);
else
    disp@handle(obj);
end


end
