%Roshan Jaiswal-Ferri
%Section - 
%Aero 302 Homework 1 FS2- 9/26/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Part 1: FS1

%Summary:
%Model the viscosity of air with altitude using the stdatm

%vars:
muRef = 1.716e-5; %pa*s
Tref = 273.15; %K
Sair = 110.4; %K

%analysis:


altVector = linspace(1,100000,10000); %Creating a vector to store altitude values

for i = 1:length(altVector) %Using flow control to store the values of T, P, and rho at different altitides in vectors
    [Ta(i), Pa(i), rhoa(i)] = stdatm_Jaiswsal_FerriRoshan(altVector(i)); %calling function and using the vector of altitudes to fill the other vectors of T, P, and rho
end
for i = 1:length(altVector)
    [Talt(i),a,Palt(i),rhoalt(i),nu(i)] = atmosisa(altVector(i));
end
%use sutherlands law to find viscosity:

stdatm_muAir = (((Ta./Tref).^1.5).*((Tref+Sair)./Ta+Sair)).*muRef;

figure
plot(stdatm_muAir,altVector)
ylabel('Alt - m*10^4')
xlabel('Viscosity - Pa*S')

%%

TaC = Ta - 273.15;

figure
plot(TaC,stdatm_muAir)
xlabel('Temp C')
ylabel('Viscosity - Pa*S')

figure; %creating a figure with three graphs of Temp, Pres, and Density vs altitude
subplot(1,3,1)
plot(Ta, altVector) %'* showes the plot points with no line. two sashes (--) is a dashed line
xlabel('Temperature (K)');
ylabel('Altitude (m*10^4)')
title('Temperature vs Altitude')
grid on;

subplot(1,3,2)
plot(Pa, altVector) 
xlabel('Pressure (kPa)');
ylabel('Altitude (m*10^4)')
title('Pressure vs altitude')
grid on;

subplot(1,3,3)
plot(rhoa, altVector) 
xlabel('Density (kg/m^3 )');
ylabel('Altitude (m*10^4)')
title('Density vs altitude')
grid on;

%Results:
%The viscosity here matches matlabs data and follows the temperature trends
%seen in each stage of our atmosphere. 


% figure; %creating a figure with three graphs of Temp, Pres, and Density vs altitude
% subplot(1,3,1)
% plot(Talt, altVector) %'* showes the plot points with no line. two sashes (--) is a dashed line
% xlabel('Temperature (K)');
% ylabel('Altitude (m*10^4)')
% title('Temperature vs Altitude')
% grid on;
% 
% subplot(1,3,2)
% plot(Palt, altVector) 
% xlabel('Pressure (kPa)');
% ylabel('Altitude (m*10^4)')
% title('Pressure vs altitude')
% grid on;
% 
% subplot(1,3,3)
% plot(rhoalt, altVector) 
% xlabel('Density (kg/m^3 )');
% ylabel('Altitude (m*10^4)')
% title('Density vs altitude')
% grid on;
