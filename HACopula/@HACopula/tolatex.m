function out = tolatex(obj, type)
%TOLATEX - returns a latex expression of cdf or pdf of the HAC in *obj*,
% i.e., the input *type* should be in {'cdf', 'pdf'}.
%
% NOTE:
% Be aware that for d > 5, getting pdf becomes extremely time
% demanding.
% In Octave, the OctSymPy package is needed (https://github.com/cbm755/octsympy).
%
%
% Copyright 2018 Jan Gorecki

if strcmp(type, 'cdf')
    f = getcdf(obj);
elseif strcmp(type, 'pdf')
    f = getpdf(obj);
else
    error('HACopula:tolatex', 'HACopula: Unsupported type. Choose one from {''cdf'', ''pdf''}.');
end

%out = latex(simplify(f)); % simplification - uncommnent to simplify the formula
%out = latex(simplify(f,'IgnoreAnalyticConstraints',true)); % simplification - uncommnent to simplify the formula
out = latex(f);

d = obj.Root.Dim; % use Root for the case that one wants to get a child/descentant copula 
for i = 1:d
    out = regexprep(out, sprintf('%s%d%s','\\mathrm{u', i, '}'), sprintf('%s%d%s','u_{', i, '}'));
end
for i = d+1:d+length(obj.Root.Forks)
    out = regexprep(out, sprintf('%s%d%s','\\mathrm{theta', i, '}'), sprintf('%s%d%s','\\theta_{', i, '}'));
end
end