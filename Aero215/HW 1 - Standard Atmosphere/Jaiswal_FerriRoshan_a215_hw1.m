%Roshan Jaiswal-Ferri
%Aero 215 HW 1: 10/9/23

%%

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Part 1

%variables for part 1: 

h = 25001;
[T, P, rho] = stdatm_Jaiswsal_FerriRoshan(h);
%[T, P, rho] = teststdatm_Jaiswsal_FerriRoshan(h);
disp(['Temp in Kelvin: ',num2str(T),', Pressure in kPa: ',num2str(P),', Density in kg/m^3: ',num2str(rho)]);


%% Part 2

altVector = linspace(1,100000,10000);

for i = 1:length(altVector)
    [T(i), P(i), rho(i)] = stdatm_Jaiswsal_FerriRoshan(altVector(i));
end

figure;
subplot(1,3,1)
plot(T, altVector) %'* showes the plot points with no line. two sashes (--) is a dashed line
xlabel('Temperature (K)');
ylabel('Altitude (x1,000km)')
title('Temperature vs Altitude')
grid on;
%grid minor; %Adds more smaller grid lines
subplot(1,3,2)
plot(P, altVector) %'* showes the plot points with no line. two sashes (--) is a dashed line
xlabel('Pressure (kPa)');
ylabel('Altitude (x1,000km)')
title('Pressure vs altitude')
grid on;
%grid minor; %Adds more smaller grid lines
subplot(1,3,3)
plot(rho, altVector) %'* showes the plot points with no line. two sashes (--) is a dashed line
xlabel('Density (kg/m^3 )');
ylabel('Altitude (x1,000km)')
title('Density vs altitude')
grid on;
%grid minor; %Adds more smaller grid lines

























