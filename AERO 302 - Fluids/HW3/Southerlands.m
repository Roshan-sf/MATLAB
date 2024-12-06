function [mu_air] = Southerlands(T)

% Define constants
mu_refAir = 1.716e-5;   % [Pa.s]
T_refAir = 273.15;      % [K]
S_air = 110.4;          % [K]

% Sutherland law formula
    term1 = (T / T_refAir)^1.5;
    term2 = (T_refAir + S_air) / (T + S_air);

    mu_air = mu_refAir * term1 * term2;

end


