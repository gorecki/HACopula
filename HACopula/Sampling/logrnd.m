function out = logrnd(alpha, n)
%LOGRND - Sampling from the logarithmic distribution
%
% out = logrnd(alpha, n) draws a sample of n observations from the
% logarithmic distribution with alpha from (0, 1). This function implements
% the algorithm "LK" from Kemp (1981).
%
% Reference:
% Kemp, A. W. (1981), ÿEfficient Generation of Logarithmically Distributed Pseudo-Random
% Variablesÿ, Journal of the Royal Statistical Society: Series C (Applied Statistics), 30(3),
% 249ÿ253.
%
%
% Copyright 2018 Jan Gorecki

out = zeros(n,1);
for i = 1:n
    h = log(1 - alpha);
    u2 = rand(1,1); 
    if (u2 > alpha)
        out(i) = 1;
    else
        u1 = rand(1,1);
        q = -expm1(u1 * h);
        if (u2 < q^2)
            out(i) = floor(1 + log(u2)/log(q));
        elseif (u2 > q)
            out(i) = 1;
        else
            out(i) = 2;
        end
    end
end

end