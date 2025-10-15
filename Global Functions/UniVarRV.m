function [r, v] = UniVarRV(r0, v0, dt, mu)
%Algorithm 3.4 (Credit Howard Curtis): Given r0, v0, find r, v at time dt later.
% Usage:
%   [r, v] = rv_from_r0v0(r0, v0, dt, mu)
% Inputs:
%   r0  - 3x1 initial position vector (km)
%   v0  - 3x1 initial velocity vector (km/s)
%   dt  - time of flight (s)
%   mu  - gravitational parameter (km^3/s^2).
% Outputs:
%   r   - 3x1 position vector at t0+dt (km)
%   v   - 3x1 velocity vector at t0+dt (km/s)

    r0 = r0(:); v0 = v0(:);
    r0n  = norm(r0);                 % |r0|
    v0n  = norm(v0);                 % |v0|
    vr0  = dot(r0, v0)/r0n;          % radial velocity component v_r0  (Alg. 3.4 Step 1b).

    % Reciprocal semimajor axis alpha = 2/|r0| - |v0|^2/mu (Alg. 3.4 Step 1c).
    alpha = 2/r0n - (v0n^2)/mu;

    % Solve universal Kepler's equation for chi (Algorithm 3.3)
    chi = kepler_U(dt, r0n, vr0, alpha, mu);

    % Lagrange coefficients f, g and derivatives fdot, gdot via universal variables (Eqs. 3.69).
    [f, g, fdot, gdot, rmag] = f_and_g(chi, dt, r0n, alpha, mu);

    % Propagate state (Eqs. 3.67–3.68): r = f r0 + g v0; v = fdot r0 + gdot v0.
    r = f.*r0 + g.*v0;
    v = fdot.*r0 + gdot.*v0;

    % Normalize any tiny numerical imaginary parts to real
    r = real(r); v = real(v);

    % --------- Nested dependencies ---------

    function chi = kepler_U(dt, r0n, vr0, alpha, mu)
        % Solve universal Kepler's equation for chi using Newton's method (Alg. 3.3)
        sqrtmu = sqrt(mu);

        % Initial guess (Battin-style; robust across conic types)
        if abs(alpha) > 1e-12
            chi = sqrtmu*abs(alpha)*dt;
        else
            % Parabolic limit; use Barker-like guess
            h = norm(cross(r0, v0));
            s = 0.5*pi*sqrtmu*dt/(r0n);
            chi = sqrtmu*dt/(r0n); % scale with time
        end

        tol = 1e-8; maxit = 50;
        for k = 1:maxit
            z  = alpha*chi^2;
            [C, S] = stumpff(z);

            % Universal Kepler equation F(chi) = 0:
            F  = (r0n*vr0/sqrtmu)*chi^2*C + (1 - alpha*r0n)*chi^3*S + r0n*chi - sqrtmu*dt;

            % Derivative dF/dchi (standard closed form):
            dF = (r0n*vr0/sqrtmu)*chi*(1 - z*S) + (1 - alpha*r0n)*chi^2*C + r0n;

            delta = F/dF;
            chi   = chi - delta;

            if abs(delta) < tol, break; end
        end
    end

    function [f, g, fdot, gdot, rmag] = f_and_g(chi, dt, r0n, alpha, mu)
        % Lagrange coefficients and their derivatives using universal variables (Eqs. 3.69)
        z = alpha*chi^2;
        [C, S] = stumpff(z);

        f = 1 - (chi^2/r0n)*C;
        g = dt - (1/sqrt(mu))*chi^3*S;

        % New radius magnitude via r = f r0 + g v0; but for coefficients we need r = |r|
        r_vec = f.*r0 + g.*v0;
        rmag  = norm(r_vec);

        fdot =  sqrt(mu)/(rmag*r0n) * (alpha*chi^3*S - chi);
        gdot = 1 - (chi^2/rmag)*C;
    end
end