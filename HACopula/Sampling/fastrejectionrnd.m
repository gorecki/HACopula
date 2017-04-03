function out = fastrejectionrnd(v0, h, c, Fm)
%FASTREJECTIONRND - Generating a sample using the fast rejection 
% algorithm.
%
% out = fastrejectionrnd(v0, h, c, Fm) generates a sample of n (where n = size(v0,1)) observations
% using the fast rejetion algorithm introduced in [Hofert, 2010].
%
% Inputs:
% v0    - A vector containing a sample from the outer distribution F0.
% h     - A tilt (see Section 4.2.2 from [Hofert, 2010] for more details).
% c     - A constant connected to \psi((c^theta + t)^(1/theta) - c) (see
%         Algorithm 4.2.9 from [Hofert, 2010] for more details).
%
% References:
% [Hofert, 2010] Hofert, M. (2010). Sampling Nested Archimedean Copulas
%     with Applications to CDO Pricing. Suedwestdeutscher Verlag fuer
%     Hochschulschriften.
%
%
% Copyright 2017 Jan Górecki

% -------------------------------------------------------------------------

n = size(v0,1);

V_sum = zeros(n,1);

for i = 1:n

%     if log_c <= 1
%         m = 1;
%     else if floor(log_c)*c(i)^(1/floor(log_c)) <= ceil(log_c)*c(i)^(1/ceil(log_c))
%             m = floor(log_c);
%         else
%             m = ceil(log_c);
%         end
%     end

    m = max(1, floor(log(c(i)))); % asymptotics

    
    V_fast = zeros(m,1);
    for j = 1:m
        %standard rejection
        cont = true;
        while cont 
            U = rand(1, 1);
            V = Fm(v0(i), m);
            W = exp(-h * V);

            cont = 1-(U <= W);
        end
        V_fast(j) = V;
    end
    
    V_sum(i) = sum(V_fast);
end

out = V_sum;