%Roshan Jaiswal-Ferri
%Aero 215 HW 1: 10/9/23

%%

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Part 1

%variables for part 1: 
%h = Given test height for function
%T = output temperature (K) 
%P = output pressure (kPa)
%rho = output density (kg/m^3)


h = 24384;
[T, P, rho] = stdatm_Jaiswsal_FerriRoshan(h); %calling the function
disp(['Temp in Kelvin: ',num2str(T),', Pressure in kPa: ',num2str(P),', Density in kg/m^3: ',num2str(rho)]);


%% Part 2

altVector = linspace(1,100000,10000); %Creating a vector to store altitude values

for i = 1:length(altVector) %Using flow control to store the values of T, P, and rho at different altitides in vectors
    [T(i), P(i), rho(i)] = stdatm_Jaiswsal_FerriRoshan(altVector(i)); %calling function and using the vector of altitudes to fill the other vectors of T, P, and rho
end

figure; %creating a figure with three graphs of Temp, Pres, and Density vs altitude
subplot(1,3,1)
plot(T, altVector) %'* showes the plot points with no line. two sashes (--) is a dashed line
xlabel('Temperature (K)');
ylabel('Altitude (x1,000km)')
title('Temperature vs Altitude')
grid on;

subplot(1,3,2)
plot(P, altVector) 
xlabel('Pressure (kPa)');
ylabel('Altitude (x1,000km)')
title('Pressure vs altitude')
grid on;

subplot(1,3,3)
plot(rho, altVector) 
xlabel('Density (kg/m^3 )');
ylabel('Altitude (x1,000km)')
title('Density vs altitude')
grid on;

