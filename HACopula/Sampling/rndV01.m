function v01 = rndV01(family0, family1, theta0, theta1, v0)
%RNDV01 - Sampling from the inverse Laplace-Stieltjes trasform of an inner
%generator.
%
% v01 = rndV01(family0, family1, theta0, theta1, v0) returns a sample of
% n(=size(v0,1)) observation from the inverse Laplace-Stieltjes
% transform of the inner generator given by family0, family1 and the
% corresponding parameters theta0, theta1.
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

TAU_NUM_UNSTABLE_A20 = 0.89;
TAU_NUM_UNSTABLE = 0.8;


familyNames = [family0 family1];

n = size(v0, 1);
switch familyNames
    case 'AA'
        v01 = zeros(n,1);
        denominator = log(1-((1-theta1)/(1-theta0)));
        for i=1:n
            v01(i) = sum(ceil(log(rand(v0(i),1))/denominator));
        end
    case 'AC'
        v = gamrnd(v0,1-theta0, n, 1);
        St = tiltedstablernd(v.^theta1, 1/theta1, ones(n,1));
        v01 = St .* v .^ theta1;
    case 'A19'
        iTheta0 = 1 - theta0;
        theta01 = theta0 * theta1;
        Fm = @(v0, m) gamrnd(iTheta0 * gamrnd(v0/m,1)/theta01, 1);
        v01 = fastrejectionrnd(v0, exp(theta1) - 1, 1./(theta0.^v0), Fm);
    case 'A20'
        th01 = (1/theta0 - 1)^theta1;
        iTheta1 = 1/theta1;
        costh1 = cos(pi/2/theta1);
        Fm = @(v0, m) gamrnd( th01 * stablernd(iTheta1, 1, (costh1/m)^theta1, 0 , 1, 1) * ...
                (gamrnd(v0/m, 1))^theta1 , 1);
        v01 = fastrejectionrnd(v0, exp(1) - 1, 1./(theta0.^v0), Fm);
        if theta2tau(family1,theta1) >= TAU_NUM_UNSTABLE_A20
            warning(['rndV0: sampling the family combination (' family0 ',' family1 ') with the child parameter (corresponding to tau = ' num2str(theta2tau(family1,theta1)) ') >= ' num2str(TAU_NUM_UNSTABLE_A20) ' may'...
                ' not produce any sample at all due to numerical instability. Use another family combination instead.']);
        end
        if theta2tau(family1,theta1) >= TAU_NUM_UNSTABLE
            warning(['rndV0: sampling the child family ' family1 ' with the parameter (corresponding to tau = ' num2str(theta2tau(family1,theta1)) ') >= ' num2str(TAU_NUM_UNSTABLE) ' may'...
                ' produce a sample that is not uniformly distributed in its univariate margins.']);
        end  
    case 'C12'
        S = stablernd(1/theta1, 1, (cos(pi/2/theta1))^theta1, 0 , n, 1);
        St = tiltedstablernd(1, theta0, v0);
        v01 = S .* St .^ theta1;
    case 'C14'
        theta0 = theta0*theta1; % apply the same process as for 'C12' but with theta0 = theta0*theta1
        S = stablernd(1/theta1, 1, (cos(pi/2/theta1))^theta1, 0 , n, 1);
        St = tiltedstablernd(1, theta0, v0);
        v01 = S .* St .^ theta1;
    case 'C19'
        cosTheta0 = cos(pi*theta0/2);
        iTheta0 = 1/theta0;
        Fm = @(v0, m) gamrnd(stablernd(theta0, 1, (cosTheta0*v0/m)^(iTheta0)/theta1, 0 , 1, 1), 1);
        v01 = fastrejectionrnd(v0, exp(theta1) - 1, exp(v0), Fm);
    case 'C20'
        if theta2tau(family1,theta1) >= TAU_NUM_UNSTABLE
            warning(['rndV0: sampling the child family ' family1 ' with the parameter (corresponding to tau = ' num2str(theta2tau(family1,theta1)) ') >= ' num2str(TAU_NUM_UNSTABLE) ' may'...
                ' produce a sample that is not uniformly distributed in its univariate margins.']);
        end
        % apply the following transformation of the parameters
        theta0 = theta0/theta1;
        theta1 = 1;
        % and then use the 'C19' approach
        cosTheta0 = cos(pi*theta0/2);
        iTheta0 = 1/theta0;
        Fm = @(v0, m) gamrnd(stablernd(theta0, 1, (cosTheta0*v0/m)^(iTheta0)/theta1, 0 , 1, 1), 1);
        v01 = fastrejectionrnd(v0, exp(theta1) - 1, exp(v0), Fm);
    case 'CC'
        v01 = tiltedstablernd(1, theta0/theta1, v0);
    case 'FF'
        % the following two constants are set according to the file
        % \copula\man\rF01FrankJoe.Rd  in the R copula package
        REJ = 1; % according to Hofert 12 (A stochastic representation ...), REJ should be from (0, infty) according to the speed of the considered computer
        APPROX = 10000;
        
        % precompute what does not change in the loops
        alpha = theta0/theta1;
        iAlpha = (theta1-theta0)/theta1;
        p0 = 1-exp(-theta0);
        p1 = 1-exp(-theta1);

        v01 = zeros(n,1);
        for i = 1:n
            if v0(i)*theta0*p0^(v0(i)-1) > REJ
                % sample via F01 for Joe
                U = rand(1);
                v01(i) = F01Joe(v0(i), theta0, theta1, APPROX);
                while U > p1^v01(i)
                    U = rand(1);
                    v01(i) = F01Joe(v0(i), theta0, theta1, APPROX);
                end
            else
                % sample as the V0-fold sum where the summands are sampled
                % via rejection with a logarithmic envelope
                for j = 1:v0(i)
                    U = rand(1);
                    %X = logrnd(1 - exp(-theta1), 1);
                    X = logrnd(min([1 - exp(-theta1) 1-eps*1e1]), 1); % NOTE: correction for theta1 > 36 
                    while (U*(X-alpha) > 1/beta(max([X eps]), iAlpha))
                        U = rand(1);
                        %X = logrnd(1 - exp(-theta1), 1);
                        X = logrnd(min([1 - exp(-theta1) 1-eps*1e1]), 1);
                    end
                    v01(i) = v01(i) + X;
                end
            end
        end
        if theta2tau(family1,theta1) >= TAU_NUM_UNSTABLE
            warning(['rndV0: sampling the child family ' family1 ' with the parameter (corresponding to tau = ' num2str(theta2tau(family1,theta1)) ') >= ' num2str(TAU_NUM_UNSTABLE) ' may'...
                ' produce a sample that is not uniformly distributed in its univariate margins.']);
        end 
    case {'GG', '1212'};
        v01 = zeros(n,1);
        alpha = theta0/theta1;
        iAlpha = 1/alpha;
        cosTheta01 = cos(pi*theta0/(2*theta1));
        for i=1:n
            v01(i) = stablernd(alpha, 1, (cosTheta01*v0(i))^iAlpha, 0, 1, 1);
        end
    case 'JJ'
        % the constant APPROX is set according to the file
        % \copula\man\rF01FrankJoe.Rd  in the R copula package
        APPROX = 10000;

        v01 = zeros(n,1);
        for i=1:n
            v01(i) = F01Joe(v0(i), theta0, theta1, APPROX);
        end
        if theta2tau(family1,theta1) >= TAU_NUM_UNSTABLE
            warning(['rndV0: sampling the child family ' family1 ' with the parameter (corresponding to tau = ' num2str(theta2tau(family1,theta1)) ') >= ' num2str(TAU_NUM_UNSTABLE) ' may'...
                ' produce a sample that is not uniformly distributed in its univariate margins.']);
        end
    case '1919'
        v01 = tiltedstablernd(exp(theta1), theta0/theta1, v0);          
    otherwise
        error(['sample_v01: class ' familyNames ' combination not supported']);
end

end

function v01 = F01Joe(v0, theta0, theta1, approx)
% sample one sample according to Joe
if v0 > approx  % approximate by stable(theta0\theta1)
    v01 = v0^(theta1/theta0) * stablernd(theta0/theta1, 1, (cos(pi*theta0/(2*theta1)))^(theta1/theta0), 0, 1, 1);
else
    v01 = sum(sibuyarnd(theta0/theta1, v0));
end
end
