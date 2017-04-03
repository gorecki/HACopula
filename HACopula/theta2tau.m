function tauMatrix = theta2tau(family, thetaMatrix)
%THETA2TAU - Computes the tau corresponding to a theta assuming *family*.
%
% tauMatrix = theta2tau(family, thetaMatrix) returns tauMatrix that corresponds
% to the given thetaMatrix through the \tau(\theta) relationship, e.g., see Table
% 2.2 in [Hofert, 2010]. Note that the function works matrix-wise.
%
% NOTE: Due to numerical unstability of the function close to the bounds of
% the parameter ranges of some of the considered families, interpolation is
% done for the cases where theta is close to the bounds (see the code below
% for details).
%
% References:
% [Górecki et al., 2016b] Górecki, J., Hofert, M., and Holeòa, M. (2016). On
%     structure, family and parameter estimation of hierarchical
%     Archimedean copulas. arXiv preprint arXiv:1611.09225.
% [Hofert, 2010] Hofert, M. (2010). Sampling Nested Archimedean Copulas
%     with Applications to CDO Pricing. Suedwestdeutscher Verlag fuer
%     Hochschulschriften.
%
%
% Copyright 2017 Jan Górecki and Martin Holeòa

tauMatrix = zeros(size(thetaMatrix));

for i = 1:size(thetaMatrix,1)
    for j = 1:size(thetaMatrix,2)
        tauMatrix(i,j) = computeKendallTau(thetaMatrix(i,j), family);
    end
end

end

%------------------------------------------------------------------------
function out = computeKendallTau(theta, family)
%returns kendall tau for parameter theta of generator from family

switch family
    case 'A'
        LOWER_BOUND_A = 1.0000e-6; % numerically unstable below the bound
        if theta >= LOWER_BOUND_A
            % standard computation
            out = 1 - 2*(theta + (1 - theta)^2 *log(1 - theta))/(3*theta^2);
        else
            % linear interpolation
            out = theta2tau('A', LOWER_BOUND_A) * theta / LOWER_BOUND_A;
        end
    case 'C'
        out = theta/(theta+2);
    case 'F'
        if theta >= 1e-4
            out = 1 - 4*(1-integral(@integrandFrank, 0, theta)/theta)/theta;
        else
            out = theta*computeKendallTau(1e-4,'F')/1e-4;
        end
    case 'G'
        out = (theta - 1)/theta;
    case 'J'
        PRECISION = 10000; % see [Hofert, 2010]
        tau = zeros(1,PRECISION);
        for k = 1:PRECISION
            tau(k) = 1/( (k*(theta*k + 2)*(theta*(k-1)+2) ));
        end
        out = max(1-4*sum(tau), 0);
    case '12'
        out = 1 - 2/(3*theta);
    case '14'
        out = 1 - 2/(1 + 2*theta);
    case '19'
        UPPER_BOUND_19 = 91;
        LOWER_BOUND_19 = 1.0e-14; % numerically unstable below the bound
        if theta > UPPER_BOUND_19
            DELTA = 1.0e-6;
            x = [UPPER_BOUND_19-DELTA UPPER_BOUND_19];
            y = [theta2tau('19', UPPER_BOUND_19 - DELTA) theta2tau('19', UPPER_BOUND_19)];
            out = min(1-eps(1), interp1(x, y, theta, 'linear', 'extrap'));
        elseif theta >= LOWER_BOUND_19
            RIGHT_BOUND = 100;
            s = warning('off', 'MATLAB:integral:MinStepSize');  % switch of the warning
            out = 1/3 + 2*theta*(1 - theta * exp(theta) * ...
                integral(@(x)exp(-x)./x, min(theta,RIGHT_BOUND), RIGHT_BOUND))/3;
            warning(s);
        else
            % linear interpolation
            out = (theta2tau('19', LOWER_BOUND_19) - 1/3) * theta / LOWER_BOUND_19 + 1/3;
        end
    case '20'
        % introduced in [Górecki et al., 2016b]
        LOWER_BOUND_20 = 1.0e-8; % numerically unstable below this bound
        if theta >= LOWER_BOUND_20
            % standard computation
            innerInteg = integral(@(t)((t.^(theta+1))./(exp(t.^(-theta)))), 0, 1);
            out = 1 - (4/theta) * (1/(theta+2) - exp(1) * innerInteg );
        else
            % linear interpolation
            out = theta2tau('20', LOWER_BOUND_20) * theta / LOWER_BOUND_20;
        end
    case '?'  
        out = theta; % an arbitrary identity tau=theta for the ? family
    otherwise
        error('theta2tau: computeKendallTau: unknown family');
end
end

function output = integrandFrank(input)

output(exp(input)>1) = input(exp(input)>1)./(exp(input(exp(input)>1))-1);

end
