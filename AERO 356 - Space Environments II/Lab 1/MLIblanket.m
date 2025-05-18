%% Lab 1D MLI Blanket

clc
clear
close all

Ekap =.72;
EDAC = .935;
Emylar = .76;

Ein= Ekap;
Eout = Ekap;
a = 2;
b = 8;
c = 5;
n = a + b+ c +2;

Eeff = (((Ein^-1)+(Eout^-1) -1)+ (((2*a)/Ekap) + ((2*b)/EDAC) + (2*c)/Emylar)-n)^-1;

disp(['Effective Emissivity ', num2str(Eeff)]);
