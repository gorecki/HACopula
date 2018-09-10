function dist = distance(obj, objToCompare, varargin)
%DISTANCE - computes a standardized Euclidian-like distance between a
% selected type of matrix (kendall, upper-tail, lower-tail)
% corresponding to obj and to objToCompare, where objToCompare
% can be either an instance of the HACopula class or a Kendall
% correlation matrix
%
% Note:
% The type of matrix can be supplied as the third argument. If it
% is not supplied, the default value is 'kendall'. Also, if
% objToCompare is a Kendall correlation matrix, do not supply
% the third argument as the 'kendall' type will be used.
%
%
% Copyright 2018 Jan Gorecki

% Inputs checking
narginchk(2,3);
d = obj.Dim;
if isa(objToCompare,'HACopula')
    if d ~= objToCompare.Dim
        error('HACopula:distance', 'HACopula::distance: Both HACopula inputs must be of the same dimension.');
    end
    type = varargin{1};
    if sum(strcmp({'kendall', 'upper-tail', 'lower-tail'}, type)) == 0
        error('HACopula:distance', 'HACopula::distance: The input *type* must be from {''kendall'', ''upper-tail'',''lower-tail''}.');
    end
elseif ismatrix(objToCompare)
    if ~((max(size(objToCompare)) == min(size(objToCompare))) && ...
            (max(size(objToCompare)) == d))
        error('HACopula:distance', 'HACopula::distance: The matrix objToCompare must be a of size d*d.');
    end
    type = 'kendall';
    if ~isempty(varargin)
        warning('HACopula:distance', 'HACopula::distance: Ignoring the third argument and setting type=''kendall''.')
    end
else
    error('HACopula:distance', 'HACopula::distance: The input objToCompare must be an instance of HACopula class or a Kendall correlation matrix.');
end

% Distance computation
dist = 0;
for i = 1:d
    for j = i+1:d
        [~, fork] = getyca(obj, i, j);
        if isa(objToCompare,'HACopula')
            [~, forkToCompare] = getyca(objToCompare, i, j);
            switch type
                case 'kendall'
                    dist = dist + (fork.Tau - forkToCompare.Tau)^2;
                case 'upper-tail'
                    dist = dist + ...
                        (gettaildependence(fork.Family, fork.Parameter, 'upper') - ...
                        gettaildependence(forkToCompare.Family, forkToCompare.Parameter, 'upper'))^2;
                case 'lower-tail'
                    dist = dist + ...
                        (gettaildependence(fork.Family, fork.Parameter, 'lower') - ...
                        gettaildependence(forkToCompare.Family, forkToCompare.Parameter, 'lower'))^2;
            end
        else % objToCompare is a matrix d*d
            % type is kendall
            dist = dist + (fork.Tau - objToCompare(i,j))^2;
        end
    end
end
% standardization
dist = sqrt(dist/nchoosek(d,2));
end