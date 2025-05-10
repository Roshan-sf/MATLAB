%% Lab 1D MLI Blanket

clc
clear
close all

%Effective Emissivity 

Ekap =.72;
EDAC = .935;
Emylar = .76;
m = 2.862; %grams

Ein= Ekap;
Eout = Ekap;
a = 2;
b = 8;
c = 5;
n = a + b+ c ;

Eeff = (((Ein^-1)+(Eout^-1) -1)+ (((2*a)/Ekap) + ((2*b)/EDAC) + (2*c)/Emylar)-n)^-1;

disp(['Effective Emissivity ', num2str(Eeff)]);

% Calultating Temperature Variation Per Unit Mass
deltaTb = 120.6 - 22;
deltaTt = 26-23.1;
tv_per_unitmass = (1/m)*(deltaTb- deltaTt);
disp(['Calultating Temperature Variation Per Unit Mass ', num2str(tv_per_unitmass)]);
