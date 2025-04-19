%% Roshan Jaiswal-Ferri
%Section - 01
%Aero 421 HW1: 4/2/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% 12.1

phi = pi/4;
theta = pi/4;
psi = pi/4;

cpi = cos(phi);
spi = sin(phi);
ct = cos(theta);
st = sin(theta);
cps = cos(psi);
spc = sin(psi);

Cbg = [ ...
    ct * cps,                      ct * spc,                          -st;
    spi * st * cps - cpi * spc,    spi * st * spc + cpi * cps,   spi * ct;
    cpi * st * cps + spi * spc,    cpi * st * spc - spi * cps,   cpi * ct
];

%disp(num2str(Cbg))

mu = 398600;
r = [0;0;7000];
rb = Cbg*r;
rbm = norm(rb);
rcross = vcross(rb);

I = [100,0,0;...
    0,120,0;...
    0,0,80];

I = I./1e6;

Tg = ((3*mu)/(rbm^5)) * rcross*I*rb;

disp('Gravity Gradient Torque (N-km): ')
disp(num2str(Tg))
disp(' ')

%% 13.1
% a) the satellite exhibits free precession with nutation due to imperfect 
% deployment, but remains stable due to spinning about the axis of max inertia.


%Assume It = Ix =Iy & Iz = Ia, where It = 100

wx = 0.1;
wy = 0.02;
wz = 0.5;

Ix = 98; % kg·m2
Iy = 102; % kg·m2
Iz = 150; %kg·m2

It = mean([Ix,Iy]);
Ia = Iz;
    
I = [It,0,0;...
    0,It,0;...
    0,0,Ia];

hx = It * wx;
hy = It * wy;
hz = Ia * wz;

ht = sqrt((hx^2) + (hy^2));
h = norm([hx,hy,hz]);

nutation = asin(ht/h);
nutation2 = acos(hz/h); % they match, yay!

processionRate = h/It;

disp('Nutation angle (rad): ')
disp(num2str(nutation))
disp(' ')
disp('Nutation angle (rad/s): ')
disp(num2str(processionRate))

%% Functions

function [out] = vcross(v)
    out = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
end