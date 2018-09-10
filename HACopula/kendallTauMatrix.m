function tau = kendallTauMatrix(data)
%KENDALLTAUMATRIX - computes Kendall's tau (if available with a n log(n) C++ implementation)
% 
% tau = kendallTauMatrix(data) - computes the matrix of Kendall's pairwise
% taus for data. If a .mex file has been compiled from fastKendallTau.cpp,
% the computation is in O(n*log(n)), otherwise in O(n^2).
%
% Input:
% data  - A (n x d)-matrix of doubles.
%
% Output:
% tau   - A (d x d)-matrix of Kendall's pairwise taus.
%
% To compile fastKendallTau.cpp (requires a C++ compiler, see, e.g.,
% www.mathworks.com/support/compilers.html), execute in the command line
% the following:
% (for MATLAB)
% mex fastKendallTau.cpp
%
% (for Octave)
% mkoctfile --mex fastKendallTau.cpp
% 
%
% Copyright 2018 Matle Kurz & Jan Gorecki

if isoctave
    if exist('fastKendallTau', 'file')==3 % check whether a compiled mex file exists
        tau = fastKendallTau(data);
    else
        tau = kendall(data);
    end
else % MATLAB
    if exist('fastKendallTau', 'file')==3 % check whether a compiled mex file exists
        tau = fastKendallTau(data);
    else
        tau = corr(data, 'type', 'kendall');
    end
end
    
end

