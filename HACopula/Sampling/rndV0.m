function v0 = rndV0(family, theta, n)
%RNDV0 - Sampling from the inverse Laplace-Stieltjes transform of an Archimedean
% generator
%
% v0 = rndV0(family, theta, n) draws a sample of n observation from the
% inverse Laplace-Stieltjes trasform of the Archimedean generator from
% family with the parameter theta.
%
% Inputs:
% family    - Family label of the generator. Possible families are: 'A', 
%             'C', 'F', 'G', 'J', '12', '14', '19' and '20'.
% theta     - The parameter of the generator.
% n         - The number of observations to be generated.
% 
% Output:
% v0        - A vector of n observations.
%
% References:
% [Hofert, 2010] Hofert, M. (2010). Sampling Nested Archimedean Copulas
%     with Applications to CDO Pricing. Suedwestdeutscher Verlag fuer
%     Hochschulschriften.
%
%
% Copyright 2017 Jan Górecki

% NOTE: for tau close to 1, the sampling procedure shows numerical problems
% several families, which is addressed by the constants below

TAU_NUM_UNSTABLE = 0.8;
TAU_NUM_UNSTABLE_F = 0.89;

switch family
    case 'A'
        v0 = ceil(log(rand(n,1))/log(theta));
    case 'C'
        v0 = gamrnd(1/theta,1,n,1);
    case 'F'
        v0 = logrnd(1 - exp(-theta), n);
        if theta2tau(family,theta) >= TAU_NUM_UNSTABLE_F
            warning(['rndV0: sampling the family ' family ' with the parameter (corresponding to tau = ' num2str(theta2tau(family,theta)) ') >= ' num2str(TAU_NUM_UNSTABLE_F) ' may'...
                ' not produce any sample at all due to numerical instability. Use another family instead.']);
        end 
    case 'G'
        v0 = stablernd(1/theta, 1, (cos(pi/(2*theta)))^theta, 0 , n, 1);
    case 'J'
        v0 = sibuyarnd(1/theta, n);
        if theta2tau(family,theta) >= TAU_NUM_UNSTABLE
            warning(['rndV0: sampling the family ' family ' with the parameter (corresponding to tau = ' num2str(theta2tau(family,theta)) ') >= ' num2str(TAU_NUM_UNSTABLE) ' may'...
                ' produce a sample that is not uniformly distributed in its univariate margins.']);
        end        
    case '12'
        v0 = stablernd(1/theta, 1, (cos(pi/(2*theta)))^theta, 0 , n, 1) .* (exprnd(1, n, 1).^theta);       
    case '14'
        v0 = stablernd(1/theta, 1, (cos(pi/(2*theta)))^theta, 0 , n, 1) .* (gamrnd(theta, 1, n, 1).^theta);
    case '19'
        v0 = gamrnd(exprnd(1, n, 1)./theta, 1/exp(theta)); 
    case '20'
        v0 = gamrnd(gamrnd(1/theta, 1, n, 1), 1/exp(1));
        if theta2tau(family,theta) >= TAU_NUM_UNSTABLE
            warning(['rndV0: sampling the family ' family ' with the parameter (corresponding to tau = ' num2str(theta2tau(family,theta)) ') >= ' num2str(TAU_NUM_UNSTABLE) ' may'...
                ' produce a sample that is not uniformly distributed in its univariate margins.']);
        end
    otherwise
        error(['rndV0: family ' family ' is not supported']);
end

