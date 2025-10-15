function [r_ECI, v_ECI, r_PERI, v_PERI] = COES2RV(ecc, theta, inc, RAAN, omega, r_p)
% Takes COEs and outputs R and V in ECI and Perifocal Frames
%   ecc = eccentricity
%   theta = true anomaly, rad
%   inc = inclination, rad
%   RAAN = right angle of ascending node, rad
%   omega = argument of perigee, rad
%   r_p = radius of perigee, km

mu = 398600; % km^3/s^2, earth gravitational constant

% find h from orbit equation
h = sqrt(mu*(1 + ecc)*r_p); % specific angular momentum, km^2/s

% find r and v relative to perifocal frame
r_PERI = h^2/mu * 1/(1+ecc*cos(theta)) * [cos(theta); sin(theta); 0]; % r in perifocal frame, km
v_PERI = mu/h * [-sin(theta); ecc+cos(theta); 0]; % v in perifocal frame, km/s

% find r and v relative to geocentric equatorial frame
[DCM_ECI_PERI,~] = DCM_PERI_ECI(omega,inc,RAAN); % call directional cosine function to transform from perifocal to ECI frame
r_ECI = DCM_ECI_PERI * r_PERI; % r in ECI frame, km
v_ECI = DCM_ECI_PERI * v_PERI; % v in ECI frame, km/s

end