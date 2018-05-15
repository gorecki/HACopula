function tau = kendallTauMatrix(data)
%kendallTauMatrix Compute Kendall's tau (if available with a n log(n) C++ implementation)

if isoctave
    tau = kendall(data(indices,:));
else % MATLAB
    if exist('fastKendallTau', 'file')==3 % check whether a compiled mex file exists
        tau = fastKendallTau(data);
    else
        tau = corr(data, 'type', 'kendall');
    end
end
    
end

