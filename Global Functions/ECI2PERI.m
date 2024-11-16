function [EtoP, PtoE] = ECI2PERI(omega,inc,RAAN)
    %ECI2PERI Earth Centered Inertial Frame to Perifocal Frame
    %   [EtoP, PtoE] = ECI2PERI(omega,inc,RAAN)
    EtoP = inv(rotz(omega))*inv(rotx(inc))*inv(rotz(RAAN));
    PtoE = inv(EtoP);
end