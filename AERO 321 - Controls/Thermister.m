%Roshan Jaiswal-Ferri
%Section - 01
%Aero 321 In Class - Lab 1: 9/25/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: Two Body Equation of Motion
filename = "ThermistorData3.csv";
filename2 = "ThermistorData5";
data = readmatrix(filename);
data2 = readmatrix(filename2);

%cleaning data
tC = data(1,1) - 1;
data(:,1) = data(:,1) - tC;
data(:,1) = data(:,1)./1000;

%cleaning data
tC = data2(1,1) - 1;
data2(:,1) = data2(:,1) - tC;
data2(:,1) = data2(:,1)./1000;

data2l = data2(1:300,:);

%Vars

r0 = 10.84; %@74 F / 23 C
r = 10000; %resister size
bVal = 3950; %Kelvin
cnst = 5/1023;

ADC = data(:,2);
dataR = resistance(ADC,r);
dataT(:,2) = RtoT(bVal,r0,dataR);

%Room Temp: 75 F (23.8889 C)

avgData = data();

%% Finding Time const

tVal = max(data2l(:,2)) * .632;

for h = 1:300
    if tVal <= data2l(h,2)
        time_constant = data2l(h,2);
        break
    end
end

disp(num2str(time_constant));


%%

figure
plot(data(:,1),data(:,2))
xlabel('Time (S)')
ylabel('Temp (C)')

figure
plot(data2l(:,1),data2l(:,2))
xlabel('Time (S)')
ylabel('Temp (C)')

%% Loop for finding avg
avg = 0;
for i = 1:11:110
    %disp(i)
    %disp(num2str(data(i,2)))
    avg = avg + data(i,2);
end

avg =  avg / 10;
disp(num2str(avg));


A = avg; %Therm temp (C)
B = 23.8889; % Ref temp (C)

percentDifference = (abs(A - B) / ((A + B) / 2)) * 100;

disp([num2str(percentDifference), '%'])

%% Finding resolution: ! = Efsr/2^m

M = 10;
Efsr = 5 - 0; %5v - 0 V

Q = (Efsr)/2^M;

%% Sensitivity Gain

T = 296.95; %in K
B = 3950;
R = 10000; % ohms

sens = -t^2/B*R;

%% Functions

function [Temp] = RtoT(B,r0,R)
    T0 = 298.15; %in Kelvin

    temp = (1/T0)+(1/B)*log(R./r0);
    Temp = 1./temp;
end


function [Resistance] = resistance(ADC,R)
    Resistance = R./((1023./ADC)-1);
end





