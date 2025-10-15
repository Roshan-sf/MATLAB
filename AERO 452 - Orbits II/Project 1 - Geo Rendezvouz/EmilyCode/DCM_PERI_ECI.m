function [DCM_ECI_PERI,DCM_PERI_ECI] = DCM_PERI_ECI(omega,inc,RAAN)
% DCM between PERI and ECI
%   Transform between the geocentric and perifocal frames by outputting a
%   directional cosine matrix (both directions for transformation).

% R1(inc)
R1_inc = [1 0 0; 0 cos(inc) sin(inc); 0 -sin(inc) cos(inc)];

% R3(omega)
R3_omega = [cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1];

% R3(RAAN)
R3_RAAN = [cos(RAAN) sin(RAAN) 0; -sin(RAAN) cos(RAAN) 0; 0 0 1];

% Generate transformation matrices
DCM_PERI_ECI = R3_omega * R1_inc * R3_RAAN; % from ECI to perifocal frame
DCM_ECI_PERI = DCM_PERI_ECI'; % from perifocal to ECI frame

end