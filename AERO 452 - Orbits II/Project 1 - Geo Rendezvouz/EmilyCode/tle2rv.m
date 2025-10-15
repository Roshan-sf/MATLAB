function [r_ECI, v_ECI, r_p, r_a,T] = tle2rv(n,ecc,Me, inc, RAAN, omega)
% This function takes values from TLE and outputs the R and V vectors
% assuming input mean notion (n) is given as [rev/day]
% Me, w, inc and RAAN all given as [deg]
%given:

muearth = 398600; %km^3/s^2
radius_earth = 6378; %km 

%find period, T
n=n*2*pi; %rad/s
T = 2*pi/n; %secs;
%semi-major axis km
a = (T*sqrt(muearth)/(2*pi))^(2/3);
%get parameter
p = a*(1-ecc^2);
%get h
h = (muearth*p)^(1/2);

% set inputs for iterative scheme
f = @(E) Me-E+ecc*sin(E); % initial function
df = @(E) -1+ecc*cos(E); % derivative of function
TOL = 10^(-8); % tolerance for iterative scheme

% find initial guess B for newton solver
if Me < pi
    B = Me + ecc/2;
elseif Me > pi
    B = Me - ecc/2;
end

% call newton solver to find E
% Newton's method
   count = 0;  % iteration counter
   while abs(f(B)) > TOL && count < 10000
       En = B - (f(B) / df(B));  % Newton's method iteration
       B = En;                     % Update Ei for the next iteration
       count = count + 1;                   % Increment the iteration counter
   end
  
   % Return the final value of eccentric anomaly
   E = B;

%solve for TA using tan function
theta = 2 * atan2(sqrt(1 + ecc) * sin(E / 2), sqrt(1 - ecc) * cos(E / 2));
%solve for R and V
matrix_1 = [cos(theta); sin(theta); 0];
matrix_2 = [-sin(theta); ecc + cos(theta); 0];
%in perifocal
r_bar = (h^2/muearth) * (1/(1 + ecc*cos(theta))) * matrix_1;
v_bar = (muearth/h) * matrix_2;
%DCM Matrix
[DCM_ECI_PERI,~] = DCM_PERI_ECI(omega,inc,RAAN);
r_ECI = DCM_ECI_PERI * r_bar;
v_ECI = DCM_ECI_PERI * v_bar;
%finding rp and ra
r_p = a*(1 - ecc);
r_a = a*(1 +ecc);

end
