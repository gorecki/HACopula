function HACPdf = getpdf(obj)
%HACPFD - a symbolic representation of the PDF.
%
% In Octave, the OctSymPy package is needed (https://github.com/cbm755/octsympy).
%
%
% Copyright 2017 Jan Górecki


leaves = obj.Leaves;
HACCdf = getcdf(obj);
%differentiate HACCdf
HACPdf = HACCdf;
for i = 1:length(leaves)
    HACPdf = diff(HACPdf, sym(sprintf('u%d',leaves(i))));
end
end