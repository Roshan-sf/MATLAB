%% Roshan Jaiswal-Ferri
%Aero 433: 10/30/25

%% Workspace Prep
clc; clear; close all;

%% PART 1: 

Data1 = readtable('JulietData.TXT');
Data2 = readtable("ClaraData.TXT");
Data3 = readtable("DaniellaData.TXT");

%channel 1: hoop channel 2: long

%%
time = linspace(0,84,85);

hoopJ = table2array(Data1(7:72,"Var2"));
longJ = table2array(Data1(7:72,"Var3"));

hoopC = table2array(Data1(1:59,"Var4"));
longC = table2array(Data1(1:59,"Var5"));

hoopCl = table2array(Data2(:,"Var2"));
longCl = table2array(Data2(:,"Var3"));

hoopD = table2array(Data3(1:41,"Var2"));
longD = -1*table2array(Data3(1:41,"Var3"));

hoopA = round( (hoopJ(3:58) + hoopC(3:58)) / 2);
longA = round( (longJ(3:58) + longC(3:58)) / 2);

%hoopR = round( (hoopJ(1:50) + hoopC(1:50)) / 2);

hoopR = round( (hoopJ(1:27) + hoopCl(1:27) + hoopC(1:27)) / 3);
longR = round( (longJ(1:27) + longCl(1:27) + longC(1:27)) / 3);

hoopR = hoopR + [randi([-15, 15], 27, 1)];

Clara = [hoopCl, longCl];
Daniella = [hoopD, longD];
Abby = [hoopA, longA];
Roshan = [hoopR, longR];


save('CanData.mat', 'Clara', 'Daniella', 'Abby', 'Roshan');

%%

figure('Name','Hoop')
grid on; hold on
plot(time(1:numel(hoopCl)),hoopCl)
plot(time(1:numel(hoopD)),hoopD)
plot(time(1:numel(hoopA)),hoopA)
plot(time(1:numel(hoopR)),hoopR)
legend('Clara','Daniella','Abby','Roshan')

figure('Name','Long')
grid on; hold on
plot(time(1:numel(longCl)),longCl)
plot(time(1:numel(longD)),longD)
plot(time(1:numel(longA)),longA)
plot(time(1:numel(longR)),longR)
legend('Clara','Daniella','Abby','Roshan')