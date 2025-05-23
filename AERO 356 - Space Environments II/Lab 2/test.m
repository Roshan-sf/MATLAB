%% Daniella Luciani 
clc
clear
close all
warning off

G1data = readtable('Data\lab 2 - space enviorments(Sheet1) (1)-1.csv');
G1data2 = readtable('Data\lab 2 - space enviorments(Sheet2)-1.csv');

G2data = readtable("Data\L2G2.csv");
G2data2 = readtable("Data\L2G2A.csv");

G3data = readtable("Data\L2G3.csv");
G3data2 = readtable("Data\L2G3A.csv");

G4data = readtable("Data\L2G4.csv");
G4data2 = readtable("Data\L2G4A.csv");

group1.data = G1data;
group1.data2 = G1data2;

group2.data = G2data;
group2.data2 = G2data2;

group3.data = G3data;
group3.data2 = G3data2;

group4.data = G4data;
group4.data2 = G4data2;

figure % Create one figure for all subplots

[Pdec1, I_adjusted, I_adjusted2] = IVCurve(group1,1,1);
[Pdec2, I_adjusted, I_adjusted2] = IVCurve(group2,2,2);
[Pdec3, I_adjusted, I_adjusted2] = IVCurve(group3,3,3);
[Pdec4, I_adjusted, I_adjusted2] = IVCurve(group4,4,4);

disp(num2str(Pdec1))
disp(num2str(Pdec2))
disp(num2str(Pdec3))
disp(num2str(Pdec4))


function [Pdec, I_adjusted, I_adjusted2] = IVCurve(dataS, groupNum, subplotIndex)
    cellA = 0.001209675; % m^2
    data = dataS.data;
    data2 = dataS.data2;
    current = data{:,1};
    voltage = data{:,2};
    temp = data{:,3};
    temp_celcius = temp + 273.15;
    
    q = 1.602e-19;       % Elementary charge (C)
    k = 1.381e-23;       % Boltzmann constant (J/K)
    I0 = 1.95e-12;       % Saturation current (A/m^2)
    
    I_adjusted = current - I0 .* (exp(q .* voltage ./ (k .* temp_celcius)) - 1);

    % Use subplot instead of figure
    subplot(2,2,subplotIndex)
    plot(I_adjusted, voltage)
    xlabel('Current (Amps)')
    ylabel('Voltage (V)')
    grid on
    hold on

    current2 = data2{:,1};
    voltage2 = data2{:,2};
    temp2 = data2{:,3};
    temp_celcius2 = temp2 + 273.15;

    I_adjusted2 = current2 - I0 .* (exp(q .* voltage2 ./ (k .* temp_celcius2)) - 1);

    plot(current2, voltage2)
    title(['Clean vs Post Arcing Cell (Group ' num2str(groupNum) ')'])
    legend('Clean Cell', 'Post Arcing')

    Power = voltage .* current;
    Power2 = voltage2 .* current2;

    Pmax = max(Power);
    Pmax2 = max(Power2);

    Eff1 = ((Pmax)/(1366 * cellA)) * 100;
    Eff2 = ((Pmax2)/(1366 * cellA)) * 100;
    Pdec = abs(Eff1 - Eff2);
end
