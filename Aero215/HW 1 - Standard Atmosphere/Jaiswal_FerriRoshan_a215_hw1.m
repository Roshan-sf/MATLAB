%Roshan Jaiswal-Ferri
%Aero 215 HW 1: 10/9/23

%%

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Part 1

%variables for part 1: 

h = 1;
[T, P, rho] = stdatm_Jaiswsal_FerriRoshan(h);
disp(['Temp in Kelvin: ',num2str(T),', Pressure in kPa: ',num2str(P),', Density in kg/m^3: ',num2str(rho)]);



























