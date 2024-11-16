%% Aero 351 Homework #3
clc
close all
clear all
%% 4.15
% Part A:
mu = 398600;
periAlt = 300; %km
rEarth = 6378; %km
rP = rEarth + periAlt;
ecc = 1.5; % Hyperbolic
inc = 35; %deg
RAAN = 130; % Deg
theta = 0; %deg
argPeri = 115; %deg

a = rP/(1-ecc);
h = (a*(1-ecc^2)*mu)^(1/2);
rPer = [(h^2/mu)*(1/(1+ecc*cosd(theta)))*cosd(theta) ; (h^2/mu)*(1/(1+ecc*cosd(theta)))*sind(theta) ; 0];
vPer = [(-mu/h)*sind(theta) ; (mu/h)*(ecc+cosd(theta)) ; 0];

% Part B
[DCMeciperi] = DCMperieci(argPeri, inc, RAAN);
rEci = DCMeciperi*rPer;
vEci = DCMeciperi*vPer;


function[DCMeciperi, DCMperieci] = DCMperieci(omega, inc, ohm)
% omega = rad2deg(omega);
% inc = rad2deg(inc);
% ohm = rad2deg(ohm);
DCMperieci = rotz(omega)'*rotx(inc)'*rotz(ohm)';
DCMeciperi = DCMperieci';
end
%% 6.25
AltPeri2 = 1270;
rPer2 =AltPeri2 +rEarth;
Vperi = 9;
theta2 = 100;
h_init = Vperi*rPer2;
ecc_init = 1/(rPer2*mu/h_init^2)-1;
ecc_fin = 100;

r100mag = (h_init^2/mu)*(1/(1+ecc_init*cosd(100)));
v100tan = (h_init/r100mag);
v100rad = (mu/h_init)*ecc_init*sind(100);
flightPath = atand(v100rad/v100tan);
v100norm = sqrt(v100rad^2+v100tan^2);

h_fin = sqrt(r100mag*mu*(1+ecc_fin*cosd(100)));
vtan_fin = h_fin/r100mag;
vrad_fin = (mu/h_fin)*ecc_fin*sind(100);
vfinnorm = sqrt(vrad_fin^2+vtan_fin^2);
flightPath2 =atand(vrad_fin/vtan_fin);
deltaV = sqrt(v100norm^2+vfinnorm^2-2*v100norm*vfinnorm*cosd(flightPath2 - flightPath));

%% 6.44
rPark = 6678;
rIncrease = 6978;
incChange = 20;
% Part A (Change r first, then plane, 3 delta v's)

vPark = sqrt(mu/rPark);
eHohmann = (rIncrease-rPark)/(rIncrease+rPark);
h_hohmann = sqrt(rPark*mu*(1+eHohmann));
vHohmann1 = h_hohmann/rPark;
vHohmann2 = h_hohmann/rIncrease;

deltaV1 = abs(vHohmann1-vPark);
vIncrease = sqrt(mu/rIncrease);
deltaV2 = abs(vIncrease-vHohmann2);



deltaVinc = 2*vIncrease*sind(incChange/2);

deltaVtotA = deltaV1+deltaV2+deltaVinc;


% Part B

deltaVincrad = sqrt((vIncrease-vHohmann2)^2 + 4*vHohmann2*vIncrease*sind(incChange/2)^2);
deltaVincrad_tot = deltaVincrad+deltaV1;

% Part C

deltaVinclow = sqrt((vHohmann1-vPark)^2 + 4*vHohmann1*vPark*sind(incChange/2)^2);
deltaVinclow_tot = deltaVinclow+deltaV2;
%Note: Please review

%% 6.47
rEarth = 6780; %km
m = 1000; %kg
r0 = [436 6083 2529]; %km
v0 = [-7.340 -0.5125 2.497]; %km/s
Isp = 300; %s
T = 10; % kN Thrust
tspan = [0 89*60] ;%s
burnTime = [89*60 89*60+120];%s
g0 = 9.807/1000;


options = odeset('RelTol',1e-8,'AbsTol',1e-8);
state = [r0 v0];

[~,stateNew] = ode45(@twobodymotion,tspan,state,options,mu);

% These are the r and v vectors at the start of the burn:
r1 = [stateNew(end,1), stateNew(end,2), stateNew(end,3)];
v1 = [stateNew(end,4), stateNew(end,5), stateNew(end,6)];

state2 = [r1 v1 m]; % state at start of burn
[~,state3] = ode45(@contThrust,burnTime,state2,options, T, Isp, mu, g0);

r2 = [state3(end,1), state3(end,2), state3(end,3)];
v2 = [state3(end,4), state3(end,5), state3(end,6)];
mf = state3(end, 7);

[a, e, nu, i, Ohm, Omega, Period, h, E, Tp,epsilon, Palt, Aalt]= RV2COES_RosuStefan(v2, r2, mu, rEarth);
disp("Highest Altitude: "+Aalt);
disp("Final Mass: "+ mf)
function tstate = contThrust(time, state, T, Isp, mu, g0)

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

function dstate = twobodymotion(time,state,muEarth) %dstate is derivitve of state

    %define vars

    x = state(1);
    y = state(2);
    z = state(3);

    dx = state(4); %vel
    dy = state(5); %vel
    dz = state(6); %vel


    r = norm([x y z]);

    ddx = -muEarth*x/r^3;
    ddy = -muEarth*y/r^3;
    ddz = -muEarth*z/r^3;
    dstate = [dx; dy; dz; ddx; ddy; ddz];
end
