function [y, nNans] = nanapprox(y, x)
%NANAPPROX - replaces NaNs and Infs by a real-valued approximations
%
% Checks for NaNs (Infs) in the input *y* and tries to replace them by an appropriate
% real value. It is assumed that y = f(x), where f is a continuous
% function. If y(i) is a NaN (Inf), it is replaced by y(j), where j is the index
% of the observation (x_{j,1}, ...,  x_{j,d}) that is the closest (in the
% Euclidian distace on [0, 1]^d) to (x_{i,1}, ..., x_{i,d}) provided y(j) is a
% real number.
%
% Inputs:
% y         - An (n * 1)-matrix with some evaluation of x
% x         - An (n * d)-matrix
%
% Outputs:
% y         - The original (input) y, where each NaN or (-)Inf is replaced
%             by an approximation.
% nNans     - The nunmer of NaN or (-)Inf encountered in the original y.
%
%
% NOTE:
% Assumes at least one real value in the input y.
%
%
% Copyright 2018 Jan Gorecki

nans = isnan(y) | isinf(y);
nNans = sum(nans);
if (nNans > 0)
    if nNans == length(y)
        warning('HACopula:nanapprox', 'nanapprox: No approximation has been done as no real value in the input y has been provided.');
    else
        iNans = find(nans);
        notNansInd = setdiff(1:length(y), iNans);
        for i = 1:size(iNans,1)
            %switch method
            %case 'euclidean'
            % compute the Euclidian distance
            dist = pdist2(x(iNans(i), :), x(not(nans),:));
            
            % find the closest x
            [~,ind] = min(dist);
            closest = notNansInd(ind(1));
            
            %substitute the NaN with the y-value of the closest x
            y(iNans(i)) = y(closest);
            %end
        end
    end
end

end