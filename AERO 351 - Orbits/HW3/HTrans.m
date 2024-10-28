%Roshan Jaiswal-Ferri
%Section - 03
%Aero 302 Lab 1 - Venturi Effect and Open Jet: 10/1/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: Calculating dV 
%for hrt check h between leo and meo
mu = 398600;
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

%% Finding Orbit time

a = (rat+rpt)/2;
P = pi*sqrt((a^3)/mu);
tP = P/60; %<3 is shorter than one hour







