function [Phi_rr, Phi_rv, Phi_vr, Phi_vv, dv0plus, dv0minus] = SolveCW(n, t, r0)
%Solves the CW equations in Matrix form
%INPUTS: SolveCW(n,t,r), or (n,t), or (n*t)
%INPUTS: SolveCW(mean motion (rad/s), elapsed time (scaler number), rvec)
%INPUTS: SolveCW(n*t, r0)
%OUTPUTS: [Phi_rr, Phi_rv, Phi_vr, Phi_vv, dv0plus, dv0minus]

    if nargin == 1
        nt = n;
    elseif nargin == 2 || nargin == 3
        nt = n*t;
    end

    s = sin(nt);
    c = cos(nt);

    Phi_rr = [ 4-3*c,           0,   0;
               6*(s-nt),        1,   0;
               0,               0,   c ];

    Phi_rv = [ (1/n)*s,         (2/n)*(1-c), 0;
               (2/n)*(c-1),     (1/n)*(4*s-3*nt), 0;
               0,               0,       (1/n)*s ];

    Phi_vr = [ 3*n*s,           0,   0;
               6*n*(c-1),       0,   0;
               0,               0,  -n*s ];

    Phi_vv = [ c,               2*s, 0;
              -2*s,        4*c-3,   0;
               0,               0,   c ];
    
    if nargin == 3
        dv0plus = (Phi_rv) \ -Phi_rr*r0;
        dv0minus = Phi_vr * r0 + Phi_vv * dv0plus;
    else
        dv0plus = NaN;
        dv0minus = NaN; 
    end

end