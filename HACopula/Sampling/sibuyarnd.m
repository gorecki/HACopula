function V = sibuyarnd(alpha, n)
%SIBUYARND - Sampling from the Sibuya distribution
% V = sibuyarnd(ALPHA, N) returns N random numbers from the Sibuya
% distribution with the parameter ALPHA.
%
% References:
% [Hofert, 2012] Hofert, M. (2012). A stochastic representation and
%     sampling algorithm for nested Archimedean copulas. Journal of Statistical
%     Computation and Simulation, 82(9):1239-1255.
%
%
% Copyright 2017 Jan Górecki

gamma_1_a = gamma(1 - alpha);
V = zeros(n,1);

for i = 1:n
    
    U = rand(1, 1);
    if U <= alpha
        V(i) = 1;
    else
        xMax = 1/eps;
        Ginv = ((1-U)*gamma_1_a)^(-1/alpha);
        fGinv = floor(Ginv);
        if (Ginv > xMax)
            V(i) = fGinv;
        elseif (1-U < 1/(fGinv*beta(fGinv, 1-alpha)))
            V(i) = ceil(Ginv);
        else
            V(i) = fGinv;
        end
    end
    
end

end