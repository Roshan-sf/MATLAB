%% Roshan Jaiswal-Ferri
%Aero 452 Homework 1: 9/24/25

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1:

global vc n

Rchaser_LVLH = [0; 40; 0]; %LVLH km
Vchaser_LVLH = [0;	0;	0]; % km/s

Rtarget_ECI = [3.4537e+04; -2.3865e+04; 3.9923e+03]; %ECI km
Vtarget_ECI = [-1.4726; -2.3526; -1.3222]; %km/s

mu = 398600;   % km^3/s^2 (Earth)

%% Chief's mean motion from osculating a
r = norm(Rtarget_ECI);
v = norm(Vtarget_ECI);
a = 1 / ( 2/r - v^2/mu );      % km, vis-viva
n = sqrt(mu/a^3);              % rad/s (works well for near-circular & short arcs)
P = 2*pi/n;                    % s, one chief orbit

%% Hold One: Football

a = 40;
syms dvx
eq1 = a == 2*(dvx/n);
sol = solve(eq1,dvx);
dvx = double(sol);

Vchaser_LVLH_postBurn = [dvx;	0;	0]; % km/s

tspan = [0, P];
state = [Rtarget_ECI; Vtarget_ECI; Rchaser_LVLH; Vchaser_LVLH_postBurn]; 
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

[~,relativeOrbits] = ode45(@relativeMotion,tspan,state,options,mu);

RC = [relativeOrbits(:,7),relativeOrbits(:,8),relativeOrbits(:,9)];
VC = [relativeOrbits(:,10),relativeOrbits(:,11),relativeOrbits(:,12)];

VCnorm = vecnorm(VC, 2, 2);

figure('Name','Relative Distance')
plot(RC(:,2),RC(:,1))
hold on
plot(0, 0, '*')
grid on
xlabel('Y Dist (km)')
ylabel('X Dist (km)')

%% Hold 2: Oscilatory

Rchaser_LVLH = [0; 1; 0]; %LVLH km
Vchaser_LVLH_postBurn = [0;	0;	0]; % km/s

t = P/4;
z = 1; %km
syms zdot0
eq2 = z == (zdot0/n)*sin(n*t);
sol = solve(eq2,zdot0);
zdot0 = double(sol);

Vchaser_LVLH_postBurn = Vchaser_LVLH_postBurn + [0; 0; zdot0];

tspan = [0, P];
state = [Rtarget_ECI; Vtarget_ECI; Rchaser_LVLH; Vchaser_LVLH_postBurn]; 
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

[~,relativeOrbits] = ode45(@relativeMotion,tspan,state,options,mu);

RC = [relativeOrbits(:,7),relativeOrbits(:,8),relativeOrbits(:,9)];

figure('Name','Oscillate')
plot3(RC(:,2),RC(:,1),RC(:,3))
hold on
plot3(0, 0, 0, '*')
grid on
xlabel('Y Dist (km)')
ylabel('X Dist (km)')

%% Hold 3: Vbar

Rchaser_LVLH = [0; 0.300; 0]; %LVLH km
Vchaser_LVLH_postBurn = [0;	0;	0]; % km/s

tspan = [0, P];
state = [Rtarget_ECI; Vtarget_ECI; Rchaser_LVLH; Vchaser_LVLH_postBurn]; 
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

[~,relativeOrbits] = ode45(@relativeMotion,tspan,state,options,mu);

RC = [relativeOrbits(:,7),relativeOrbits(:,8),relativeOrbits(:,9)];
VC = [relativeOrbits(:,10),relativeOrbits(:,11),relativeOrbits(:,12)];

VCnorm = vecnorm(VC, 2, 2);

figure('Name','Relative Distance')
plot(RC(:,2),RC(:,1))
hold on
plot(0, 0, '*')
grid on
xlabel('Y Dist (km)')
ylabel('X Dist (km)')

%% 2nd Football:

a = 0.02;
syms dvx
eq1 = a == 2*(dvx/n);
sol = solve(eq1,dvx);
dvx = double(sol);

Rchaser_LVLH = [0; 0.02; 0]; %LVLH km
Vchaser_LVLH_postBurn = [dvx;	0;	0]; % km/s

tspan = [0, P];
state = [Rtarget_ECI; Vtarget_ECI; Rchaser_LVLH; Vchaser_LVLH_postBurn]; 
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

[~,relativeOrbits] = ode45(@relativeMotion,tspan,state,options,mu);

RC = [relativeOrbits(:,7),relativeOrbits(:,8),relativeOrbits(:,9)];
VC = [relativeOrbits(:,10),relativeOrbits(:,11),relativeOrbits(:,12)];

VCnorm = vecnorm(VC, 2, 2);

figure('Name','Relative Distance')
plot(RC(:,2),RC(:,1))
hold on
plot(0, 0, '*')
grid on
xlabel('Y Dist (km)')
ylabel('X Dist (km)')

%% Vbar Approach

vc = 40/P; %km/s
Rchaser_LVLH = [0; 40; 0]; %LVLH km
Vchaser_LVLH_postBurn = [0;	vc;	0]; % km/s

tspan = [0, P];
state = [Rtarget_ECI; Vtarget_ECI; Rchaser_LVLH; Vchaser_LVLH_postBurn]; 
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

[~,relativeOrbits] = ode45(@vbarapp,tspan,state,options,mu);

RC = [relativeOrbits(:,7),relativeOrbits(:,8),relativeOrbits(:,9)];
VC = [relativeOrbits(:,10),relativeOrbits(:,11),relativeOrbits(:,12)];

VCnorm = vecnorm(VC, 2, 2);

figure('Name','Relative Distance')
plot(RC(:,2),RC(:,1))
hold on
plot(0, 0, '*')
grid on
xlabel('Y Dist (km)')
ylabel('X Dist (km)')


%% Functions

function dstate = vbarapp(time,state,mu) 
%Use Column vectors!
%INPUTS: first 6 rows: [x y z dx dy dz...] in ECI, target properties
%CONTD: second 6 rows: [...x y z dx dy dz] in LVLH, relative to target
%OUTPUT: follows same convention
    
    global vc n
    %unpack for clarity (t for target c for chaser):
    tx0 = state(1); %pos
    ty0 = state(2);
    tz0 = state(3);
    tdx0 = state(4); %vel
    tdy0 = state(5);
    tdz0 = state(6);

    cx0 = state(7); %pos
    cy0 = state(8);
    cz0 = state(9);
    cdx0 = state(10); %vel
    cdy0 = state(11);
    cdz0 = state(12);
    
    rvect = [tx0 ty0 tz0]; %r and v vectors for chaser from chief
    vvect = [tdx0 tdy0 tdz0];

    rt = norm([tx0 ty0 tz0]); %r vector magnitudes
    rc = rt; %norm([cx0 cy0 cz0]);
    hc = norm(cross(rvect,vvect));

    %target
    tddx = -mu*tx0/rt^3;
    tddy = -mu*ty0/rt^3;
    tddz = -mu*tz0/rt^3;

    dstate_t = [tdx0; tdy0; tdz0; tddx; tddy; tddz];

    %chaser
    cddx0 = ((2*mu/rc^3)+(hc^2/rc^4))*cx0 - 2*(dot(vvect,rvect))*(hc/rc^4)*cy0+((2*hc)/(rc^2))*cdy0;
    cddx = cddx0 + (-2*n*vc);
    cddy = ((-mu/rc^3)+(hc^2/rc^4))*cy0 + 2*(dot(vvect,rvect))*(hc/rc^4)*cx0-2*(hc/rc^2)*cdx0;
    cddz = -(mu/rc^3)*cz0;

    dstate_c = [cdx0; cdy0; cdz0; cddx; cddy; cddz];

    dstate = [dstate_t; dstate_c];

end

