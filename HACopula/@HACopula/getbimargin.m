function biMargin = getbimargin(obj, i, j)
%GETBIMARGIN - returns a bivariate margin
%
% biMargin = getbimargin(obj, i, j) returns the bivariate Archimedean
% copula corresponding to the (i, j)-th bivariate margin of HAC obj. Note
% that biMargin is an instance of HACopula.
%
% Copyright 2018 Jan Gorecki

[~, fork] = getyca(obj, i, j);
biMargin = HACopula({ {fork.Family, fork.Parameter}, 1, 2});

end