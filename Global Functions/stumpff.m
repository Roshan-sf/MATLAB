function [C, S] = stumpff(z)
    %STUMPFF: Given Z will calculate new C(z) and S(z)
    %   [C, S] = stumpff(z)

    if z > 0
        S = (sqrt(z)-sin(sqrt(z)))/((sqrt(z))^3);
        C = (1-cos(sqrt(z)))/z;
    elseif z < 0
        S = (sinh(sqrt(-z))-sqrt(-z))/((sqrt(-z))^-3);
        C = (cosh(sqrt(-z))-1)/-z;
    elseif z == 0
        S = 1/6;
        C = 1/2;
    else
        error('stumpff broke? (not a number?)')
    end
end