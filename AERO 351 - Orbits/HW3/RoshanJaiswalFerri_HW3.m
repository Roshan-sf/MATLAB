%Roshan Jaiswal-Ferri
%Section - 01 
%Aero 351: Space Debris Removal - 11/13/24

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

R = [(((h^2)/mu)*(1/(1+ecc*cosd(theta))))*cosd(theta);...
    (((h^2)/mu)*(1/(1+ecc*cosd(theta))))*sind(theta);0];
V = [(mu/h)*-sind(theta);(mu/h)*(ecc+cosd(theta));0];

[~,PtoE] = ECI2PERI(w,inc,RAAN);

R_eci = PtoE*R;
V_eci = PtoE*V;

%% PART 2: 5.6

R = [5644, 2830, 4170];
R1 = [-2240, 7320, 4980];

dt = 20*60;
tm =1;
[V0,V] = AERO351lambertUV(R,R1,tm,dt);

%% Functions:

function [EtoP, PtoE] = ECI2PERI(omega,inc,RAAN)
    EtoP = inv(rotz(omega))*inv(rotx(inc))*inv(rotz(RAAN));
    PtoE = inv(EtoP);
end

function [V0,V] = AERO351lambertUV(R0,R,tm,dt)
    mu=398600;
    tol = 1e-8;
    r0 = norm(R0);
    r = norm(R);
    
    V0=[0 0 0];
    V=[0 0 0];
    cosdnu= dot(R0,R)/(r0*r);
    A = tm*sqrt(r0*r*(1.0+cosdnu));
    
    psiold = 0;
    c2 = 0.5;
    c3 = 1/6;
    psiupper= 4.0*pi^2;
    psilower= -4.0*pi^2;
    count = 0;
    dtn = -10.0; %took from Vallado

    while abs(dtn-dt) >= tol

        y= r0 + r + ( A*(psiold*c3-1)/sqrt(c2) );
        Xn= sqrt(y/c2);
        dtn = ((Xn^3)*c3+A*sqrt(y))/sqrt(mu);

        if ( dtn < dt )
            psilower= psiold;
        end
        if ( dtn > dt )
            psiupper= psiold;
        end

        psinew= 0.5*(psiupper+psilower);
        
        [c2,c3] = findc2c3(psinew);
        psiold = psinew;
        %pause;
        
        count = count + 1;
        
        f = 1- y/r0;
        gdot= 1 - y/r;
        g = 1/(A*sqrt(y/mu));
        
        for i= 1 : 3
            V0(i)= (R(i)-f*R0(i))*g;
            V(i) = (gdot*R(i)-R0(i))*g;
        end
    end
end