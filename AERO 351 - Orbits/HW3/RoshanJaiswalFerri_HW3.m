%Roshan Jaiswal-Ferri
%Section - 01 
%Aero 351 Homework 2: 11/13/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: 4.15

mu = 398600;
Rearth = 6378;

ecc = 1.5;
Rp = Rearth + 300;
inc = 35;
RAAN = 130;
w = 115; %arg of peri
theta = 0; %at perigee

a = Rp/(1-ecc); %km
h = (mu*(a*(1-ecc^2)))^(1/2); %km2/s

R1 = [(((h^2)/mu)*(1/(1+ecc*cosd(theta))))*cosd(theta);...
    (((h^2)/mu)*(1/(1+ecc*cosd(theta))))*sind(theta);0];
V = [(mu/h)*-sind(theta);(mu/h)*(ecc+cosd(theta));0];

[~,PtoE] = ECI2PERI(w,inc,RAAN);

R_eci = PtoE*R1;
V_eci = PtoE*V;

disp('Problem 4.15:')
disp('Perifocal Frame Position:')
disp(['Mag (km): ', num2str(norm(R1))])
disp(['Vector (km): ', num2str(R1')])
disp('Perifocal Frame Velocity:')
disp(['Mag (km/s): ', num2str(norm(V))])
disp(['Vector (km/s): ', num2str(V')])
disp('Geocentric Frame Position:')
disp(['Mag (km): ', num2str(norm(R_eci))])
disp(['Vector (km): ', num2str(R_eci')])
disp('Geocentric Frame Velocity:')
disp(['Mag (km/s): ', num2str(norm(V_eci))])
disp(['Vector (km/s): ', num2str(V_eci')])
disp(' ')

%% PART 2: 5.6

R1 = [5644, 2830, 4170];
R2 = [-2240, 7320, 4980];

%parameters for lambert w/ uv function:
tol = 1e-8;
dt = 20*60;
tm = 1;

[V1,V2] = lambUVBi(R1,R2,dt,tm,mu,tol);

disp('Problem 5.6:')
disp('Starting Velocity:')
disp(['Mag (km/s): ', num2str(norm(V1))])
disp(['Vector (km/s): ', num2str(V1)])
disp('Ending Velocity:')
disp(['Mag (km/s): ', num2str(norm(V2))])
disp(['Vector (km/s): ', num2str(V2)])
disp(' ')

%% PART 3: 6.8

%Calculating dV 
mu = 398600; %for hrt check: h between leo and meo
rpt = 6378+300; %km
rat = 6378+3000; %km

V1 = sqrt(mu/rpt); %<3 is between 7 and 8!
V3 = sqrt(mu/rat);

ecc = (rat-rpt)/(rat+rpt); %<3 should be between 0 and 1!
ht = sqrt(rpt*mu*(1+ecc*cos(0))); %<3 check: h is between leo and meo

Vpt = ht/rpt; %<3 is faster than apogee
Vat = ht/rat;

dVp = Vpt - V1; 
dVa = V3 - Vat;
dV = dVa + dVp;

%Finding Orbit time
a = (rat+rpt)/2;
P = pi*sqrt((a^3)/mu);
tP = P/60; %<3 is shorter than one hour

disp('Problem 6.8:')
disp(['Required delta V (km/s): ', num2str(dV)]);
disp(['Transfer Orbit Period (min): ', num2str(tP)]);
disp(' ')

%% Reset Workspace Vars
clear all;      %Clears Workspace
mu = 398600;
Rearth = 6378;

%% PART 4: 6.23

%Orbit 1:
Rp1 = 8100; %km
Ra1 = 18900; %km
theta1B = 45; %deg
theta1C = 150; %deg
ecc1 = (Ra1-Rp1)/(Ra1+Rp1);
a1 = (Rp1+Ra1)/2;
h1 = sqrt(mu*a1*(1-ecc1^2));
T1 = 2*pi*sqrt((a1^3)/mu);

Rb1 = ((h1^2)/mu)/(1+ecc1*cosd(theta1B)); %current position of s/c B
Vb1t = (mu/h1)*(1+ecc1*cosd(theta1B)); %V tangential
Vb1r = (mu/h1)*ecc1*sind(theta1B); %V radial
Vb1 = sqrt((Vb1t^2)+(Vb1r^2));
FPAb1 = atand(Vb1r/Vb1t); %flight path angle in deg

%For s/c B
E = 2*atan(sqrt((1-ecc1)/(1+ecc1))*tand(theta1B/2)); %eccentric anomaly
Me = E-ecc1*sin(E); %Mean anomaly
ts = (Me*T1)/(2*pi); %time since periapsis

%For s/c C
E2 = 2*atan(sqrt((1-ecc1)/(1+ecc1))*tand(theta1C/2)); %eccentric anomaly
Me2 = E2-ecc1*sin(E2); %Mean anomaly
ts2 = (Me2*T1)/(2*pi); %time since periapsis

%Orbit 2:
dt = T1-(ts2-ts);
a2 = (dt^(2/3)*mu^(1/3))/(2^(2/3)*pi^(2/3));
Rp2 = Rb1; %treating position of s/c as periapsis
Ra2 = 2*a2-Rp2; 
ecc2 = (Ra2-Rp2)/(Ra2+Rp2);
h2 = sqrt(mu*a2*(1-ecc2^2));
Vb2 = h2/Rp2;

%DeltaV EQ:
%2nd Flight path angle is 0:
dV = 2*sqrt((Vb1^2)+(Vb2^2)-2*Vb1*Vb2*cosd(0-FPAb1));
% It looks like the answer is twice the dV equation, probably once for
% speeding up and once for slowing down

disp('Problem 6.23:')
disp(['Total delta V (km/s): ', num2str(dV)])
disp(' ')

%% Reset Workspace Vars
clear all;      %Clears Workspace
mu = 398600;
Rearth = 6378;

%% PART 5: 6.25

disp('done');

%% Functions:

function [EtoP, PtoE] = ECI2PERI(omega,inc,RAAN)
    EtoP = inv(rotz(omega))*inv(rotx(inc))*inv(rotz(RAAN));
    PtoE = inv(EtoP);
end

function [V1,V2] = lambUVBi(R1,R2,dtime,Tm,mu,tol)
    R1n = norm(R1);
    R2n = norm(R2);
    Z = 0;
    C = 1/2;
    S = 1/6;
    Zu = 4*pi^2; %upper bound
    Zl = -4*pi^2; %lower bound
    dtl = 1; %change in time of loop (random guess)
    A = Tm*sqrt(R1n*R2n*(1+(dot(R1,R2)/(R1n*R2n))));
    
    while abs(dtl-dtime) >= tol
        Y = R1n+R2n+(A*((Z*S-1)/sqrt(C)));        
        UV = sqrt(Y/C);
        dtl = (((UV^3)*S)/sqrt(mu))+((A*sqrt(Y))/sqrt(mu));

        if dtl < dtime
            Zl = Z; %reset zlower
        elseif dtl > dtime
            Zu = Z; %reset zupper
        end

        Z = 0.5*(Zu+Zl); %update z 
        [C,S] = stumpff(Z); %update stumpff c(z) s(z)

        f = 1-(((UV^2)/R1n)*C);
        g = (dtl - (((UV^3)/sqrt(mu)) * S))^-1;
        % fd = (sqrt(mu)/(R1n*R2n))*UV*((Z*S)-1);
        gd = 1-(((UV^2)/R2n)*C);

        %g = (1/mu)*(((Y/C)^(3/2))*S+(A*sqrt(Y)))-(1/mu)*((Y/C)^(3/2))*S;
        %g = A*sqrt(Y/mu);

        for i = 1:3
            V1(i) = (R2(i)-f*R1(i))*g;
            V2(i) = (gd*R2(i)-R1(i))*g;
        end     
    end
end

function [C, S] = stumpff(z)
    if z > 0
        S = (sqrt(z)-sin(sqrt(z)))/((sqrt(z))^3);
        C = (1-cos(sqrt(z)))/z;
    elseif z < 0
        S = (sinh(sqrt(-z))-sqrt(-z))/((sqrt(-z))^-3);
        C = (cosh(sqrt(-z))-1)/-z;
    elseif z == 0
        S = 1/6;
        C = 1/2;
    else
        error('stumpff broke? (not a number?)')
    end
end