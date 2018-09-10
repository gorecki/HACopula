function out = tiltedstablernd(h, alpha, v0)
%TILTEDSTABLERND - Sampling from the tilted stable distribution.
%
% out = tiltedstablernd(h, alpha, v0) returns a sample of
% n(=size(v0,1)) observation from the tilted stable distribution with the
% tilt h, the parameter alpha and using the outer distribution v0. The
% function implements the algorithm introduced in [Hofert, 2010].
%
% References:
% [Hofert, 2010] Hofert, M. (2010). Sampling Nested Archimedean Copulas
%     with Applications to CDO Pricing. Suedwestdeutscher Verlag fuer
%     Hochschulschriften.
%
%
% Copyright 2018 Jan Gorecki

if alpha == 1
    out = v0;
    return;
end

n = size(v0,1);

% check if h is of the same size as v0
if size(h,1) ~= n  % if not
    if size(h,1) ~=1
        error('HACopula:BadInputs', 'tiltedstablernd: parameter h must be either scalar or of the same size as v0.');
    end
    % make it of the same size
    h = ones(n,1)*h;
end
    
V_sum = zeros(n,1);

for i = 1:n
    m = max(1, round(v0(i)*h(i)^alpha));
    V_fast = zeros(m,1);
    c = (v0(i)/m)^(1/alpha); % from Marius
    for j = 1:m
        %standard rejection
        cont = true;
        
        while cont 
            U = rand(1, 1);
            % new V
            % sample V ~ S(alpha, 1, (cos(alpha*pi/2)*V_0/m)^{1/alpha},
            %   	       (V_0/m)*I_{\alpha = 1}; 1) with
            % Laplace-Stieltjes transform exp(-(V_0/m)*t^alpha)
            V = c * stablernd(alpha, 1, (cos(alpha*pi/2))^(1/alpha), 0 , 1, 1);
            W = exp(-h(i) * V);
            cont = U > W;
        end
        V_fast(j) = V;
    end
    
    V_sum(i) = sum(V_fast);
end


out = V_sum;

end