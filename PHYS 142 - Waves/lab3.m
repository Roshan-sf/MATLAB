%% Roshan Jaiswal-Ferri
%Section - 45
%PHYS 142

%% Workspace Prep

warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Plotting for part c

mf = [1,	135; % mode vs freq
      2,	277;
      3,	417;
      4,	554;
      5,	695;
      6,	787;];

x = mf(:,1);
f = mf(:,2);

p = polyfit(x,f,1);
bf = polyval(p,x);

figure('Name','Linearized Form')
plot(x,f,'*')
hold on
grid on
plot(x,bf)
% plot(0,0,'r.')
% plot(0,0,'r*')
% plot(0,0,'ro')
xlabel('Mode #')
ylabel('Freq (Hz)')
title('Mode vs Frequency')
legend('Data Points','Line of Best Fit (y = 132.88x + 12.4) ',...
    'Origin',Location='best')

syms v
eq = 277 == 2*(v/2.44);

soln = solve(eq,v);



