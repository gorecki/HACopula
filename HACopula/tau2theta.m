function thetaMatrix = tau2theta(family, tauMatrix)
%TAU2THETA - Computes the theta corresponding to a tau assuming *family*. 
%
% thetaMatrix = tau2theta(family, tauMatrix) returns thetaMatrix that corresponds
% to the given tauMatrix through the \tau(\theta) relationship, e.g., see Table
% 2.2 in [Hofert, 2010]. Note that the function works matrix-wise.
%
% NOTE: 
% Due to numerical instability of the function close to the bounds of
% [0, 1] for some of the considered families, interpolation is done for the
% cases where tau is close to the bounds of [0, 1] (see the code below for
% details).
% 
% TODO:
% In Octave, for the family F, does not work for tau > 0.9999, and for the
% family 20, does not work for tau > 0.98. An approximation (in theta2tau)
% is needed.
%
% References:
% [Hofert, 2010] Hofert, M. (2010). Sampling Nested Archimedean Copulas
%     with Applications to CDO Pricing. Suedwestdeutscher Verlag fuer
%     Hochschulschriften.
%
%
% Copyright 2018 Jan Gorecki

thetaMatrix = zeros(size(tauMatrix));

for i = 1:size(tauMatrix,1)
    for j = 1:size(tauMatrix,2)
        thetaMatrix(i,j) = computekendalltauinv(tauMatrix(i,j), family);
    end
end

end


% -------------------------------------------------------------------------

function out = computekendalltauinv(tau, family)
% Given *tau*, returns the its theta conterpart in the theta-tau
% relationship for given *family*.
% If close to the bound of \tau_{(a)}(\Theta_a), the values are linearly
% interpolated.

switch  family
    case 'A'
        UPPER_BOUND_A = 1/3 - 0.01;
        if tau >= 1/3
            theta = NaN;
        elseif tau == 0
            theta = 0;
        elseif tau <= UPPER_BOUND_A
            fce = @(t) (tau - theta2tau('A', t));
            theta = fzero(fce, [0 1-eps(1)]);
        else
            % linear extrapolation with the upper limit theta = 1-eps(1)
            DELTA = 1.0e-6;
            x = [UPPER_BOUND_A-DELTA UPPER_BOUND_A];
            y = [tau2theta('A', UPPER_BOUND_A - DELTA) tau2theta('A', UPPER_BOUND_A)];
            theta = min(interp1(x, y, tau, 'linear', 'extrap'), 1 - eps(1));
        end
    case 'C'
        theta = 2*tau/(1- tau);
    case 'F'
        %tau
        if (tau == 1)
            theta = Inf;
        else
            options = optimset('Display','off');
            fce = @(t)(tau - theta2tau(family, t));
            if (tau > 0)
                for i = 1:10
                    try
                    theta = fzero(fce, exp(i), options);
                    if ~isnan(theta) break; end
                    catch
                        % do nothing
                    end
                end
            elseif (tau < 0)
                for i = 1:10
                    theta = fzero(fce, -exp(i), options);
                    if ~isnan(theta) break; end
                end
            else
                theta = 0;
            end
        end
    case 'G'
        theta = 1/(1-tau);        
    case 'J'
        LOWER_BOUND_J = 1.0e-7; % relates to PRECISION in theta2tau.m
        UPPER_BOUND_J = 0.97;
        MAX_THETA = 1.0e+16;
        if tau < LOWER_BOUND_J
            % linear interpolation
            theta = (tau2theta('J', LOWER_BOUND_J) - 1) * tau / LOWER_BOUND_J + 1;
        elseif tau > UPPER_BOUND_J
            DELTA = 1.0e-6;
            x = [UPPER_BOUND_J-DELTA UPPER_BOUND_J];
            y = [tau2theta('J', UPPER_BOUND_J - DELTA) tau2theta('J', UPPER_BOUND_J)];
            theta = interp1(x, y, tau, 'linear', 'extrap');
        else
            % standard computation
            fce = @(t)(tau - theta2tau(family, t));
            theta = fzero(fce, [0 MAX_THETA]);
        end
    case '12'
        theta = 2/(3*(1-tau));
    case '14'
        theta = 1/(1-tau) - 0.5;
    case '19'
        UPPER_BOUND_19 = 1-1.0e-2;   
        MAX_THETA = 91;
        if tau <= 1/3
            theta = NaN;
        elseif tau > UPPER_BOUND_19
            DELTA = 1.0e-6;
            x = [UPPER_BOUND_19-DELTA UPPER_BOUND_19];
            y = [tau2theta('19', UPPER_BOUND_19 - DELTA) tau2theta('19', UPPER_BOUND_19)];
            theta = interp1(x, y, tau, 'linear', 'extrap');
        else
            fce = @(t)(tau - theta2tau('19', t));
            theta = fzero(fce, [eps(0) MAX_THETA]);
        end
    case '20'
        fce = @(t)(tau - theta2tau('20', t));
        theta = fzero(fce, tau);
    case '?'  
        theta = tau;
    otherwise
        error('HACopula:BadInputs', 'tau2theta: computekendalltauinv: unknown family');
end
out = theta;
end
