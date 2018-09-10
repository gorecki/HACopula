function HACCdf = getcdf(obj)
%HACCDF - a symbolic representation of the CDF
%
% In Octave, the OctSymPy package is needed (https://github.com/cbm755/octsympy).
%
%
% Copyright 2018 Jan Gorecki

try
    syms t;
    syms theta;
catch err
    if isoctave
        error('HACopula:OctaveToolboxNotInstalled', ['getcdf: Symbolic toolbox is needed. Please install it ' ...
            '(https://github.com/cbm755/octsympy) and load it (pkg load symbolic).']);
    else % MATLAB
        rethrow(err);
    end
end

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