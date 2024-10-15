%Roshan Jaiswal-Ferri
%Section - 
%Aero 302 Homework 1 FF1 - 9/26/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: FF1
%Summary:
%Model the viscosity using the sutherlands equation of air over the range
%of -50 to 500 degrees C


%Analysis:

%Changing viscosity due to temperature — Sutherland's Law

Temp = -50:1:500; %in C
TempK = Temp + 273.15; %in K

%muAir
muRef = 1.716e-5; %pa*s
Tref = 273.15; %K
Sair = 110.4; %K

%Sutherland's Law:

muAir = (((TempK./Tref).^1.5).*((Tref+Sair)./TempK+Sair)).*muRef;

figure
plot(Temp,muAir)
xlabel('Temp (C)')
ylabel('Viscosity (Pa*S)')

%Results:

%The fluid properties of viscocity vary between micro and macroscopic
%perspectives. On a macroscopic level viscosity can be viewed as internal
%forces that slow down a fluid. On a microscopic level viscosity is seen as
%the exchange of momentum between fluid molecules, generally slowing them
%down. When the temperature rises the molecules inside the fluid become
%excited - introducing more energy, therefore increasing the speed at which
%the fluid moves.

%From a Microscopic Perspective viscosity is dependant on pressure, the
%higher the pressure the closer molecules are together resulting in more
%collisions. On the other hand, from a macroscopic perspective, Viscocity
%is largely not dependant on pressure.

%Having a low temperature for viscosity in a wind tunnel would be important
%for testing accuracy because air at high altitudes is also very cold, so
%by lowering the temperature in the wind tunnel you are matching the
%viscosity the air would have during a real world test. 

% Yes, my plots agree with the standard atmosphere model for temperature
% pressure and viscosity. A few assumptions made by the standard atmosphere
% model and my model is 1) Hydrostatic Equilibrium 2)ideal hgas law
% 3)constant gravity 4) No atmospheric variation brought on by local
% events.

% You cannot apply exactly the same formulations to the event of a 
% de-orbiting space craft because the temperature and pressure
% magnitude and difference will be very large on an object that is 
% re-entering the atmosphere. This means
% you can no longer assume equilibrium or the ideal gas law - meaning
% behavior will be different and cannot be accurately modeled with these
% equations. What this would be useful for is predicting what the general
% atmosphere would look like and what to expect before re-entry.
















    
