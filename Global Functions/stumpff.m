function [C, S] = stumpff(z)
    % STUMPFF  Computes the Stumpff functions C(z) and S(z)
    %
    %   [C, S] = stumpff(z)
    %
    %   Inputs:
    %       z : scalar or array input
    %
    %   Outputs:
    %       C : Stumpff C(z)
    %       S : Stumpff S(z)

    if z > 1e-8
        s = sqrt(z);
        C = (1 - cos(s))/z;
        S = (s - sin(s))/(s^3);

    elseif z < -1e-8
        s = sqrt(-z);
        C = (cosh(s) - 1)/(-z);
        S = (sinh(s) - s)/(s^3);

    else
        % Series expansion around z = 0
        C = 1/2 - z/24 + z^2/720 - z^3/40320;
        S = 1/6 - z/120 + z^2/5040 - z^3/362880;
    end

end