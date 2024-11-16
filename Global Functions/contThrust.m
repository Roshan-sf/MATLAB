function [tstate] = contThrust(time, state, T, Isp, mu, g0)
    %CONTTHRUST: Propogates a new orbit trajectory during thrusting phase
    %   [tstate] = contThrust(time, state, T, Isp, mu, g0)
    %
    %   Make sure g0 (usually 9.807 m/s^2) is in km/s^2

    x = state(1);
    y = state(2);
    z = state(3);
    dx = state(4);
    dy = state(5);
    dz = state(6);
    m = state(7);
    rMag = norm([x y z]);
    vMag = norm([dx dy dz]);
    ddx = -mu*x/rMag^3 +(T/m)*dx/vMag;
    ddy = -mu*y/rMag^3 +(T/m)*dy/vMag;
    ddz = -mu*z/rMag^3 +(T/m)*dz/vMag;
    mdot = -T/(g0*Isp);
    

    tstate = [dx;dy;dz;ddx;ddy;ddz;mdot];

end