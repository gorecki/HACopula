function depMatrix = getdependencematrix(obj, type)
%GETDEPENDENCEMATRIX - returns a [d x d] matrix of a selected dependence
% coefficients (according to type from {kendall, lower-tail and
% upper-tail}). If type=tails, then tail coefficients are
% condensed into one matrix and the upper-tail coefficients are
% placed above the diagonal, whereas the lower-tail
% coefficients below it.
%
% NOTE:
% The main diagonal always contains ones as all the considered
% coefficients are 1 for two same random variables.
%
%
% Copyright 2017 Jan Górecki

if sum(strcmp({'kendall', 'upper-tail', 'lower-tail', 'tails'}, type)) == 0
    error('HACopula::getdependencematrix: The input *type* must be from {''kendall'', ''upper-tail'', ''lower-tail'', ''tails''}.')
end
d = getdimension(obj);
depMatrix = ones(d, d);
for i = 1:d
    for j = i+1:d
        [~, fork] = getyca(obj, i, j);
        switch type
            case 'kendall'
                depMatrix(i, j) = fork.Tau;
                depMatrix(j, i) = depMatrix(i, j);
            case 'upper-tail'
                depMatrix(i, j) = gettaildependence(fork.Family, fork.Parameter, 'upper');
                depMatrix(j, i) = depMatrix(i, j);
            case 'lower-tail'
                depMatrix(i, j) = gettaildependence(fork.Family, fork.Parameter, 'lower');
                depMatrix(j, i) = depMatrix(i, j);
            case 'tails'
                depMatrix(i, j) = gettaildependence(fork.Family, fork.Parameter, 'upper');
                depMatrix(j, i) = gettaildependence(fork.Family, fork.Parameter, 'lower');
        end
    end
end
end