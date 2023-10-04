function [T, P, rho] = stdatm_Jaiswsal_FerriRoshan(h)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%Variables:
%T1 is temp initial in K
%H1 is height initial in m
%a is lapse rate in gradient layers in K/m


    P1 = 101.325;
    rho1 = 1.2250;
    
    if h < 11000 && h > 0
        T1 = 288.16;
        H1 = 0;
        a = -0.0065;
        Gradient = 1;
    elseif h >= 11000 && h < 25000
        T1 = 216.66;
        H1 = 11000;
        Gradient = 0;
    elseif h >= 25000 && h < 47000
        T1 = 216.66;
        H1 = 25000;
        a = 0.003;
        Gradient = 1;
    elseif h >= 47000 && h < 53000
        T1 = 282.66;
        H1 = 47000;
        Gradient = 0;
    elseif h >= 53000 && h < 79000
        T1 = 282.66;
        H1 = 53000;
        a = -0.0045;
        Gradient = 1;
    elseif h >= 79000 && h < 90000
        T1 = 165.66;
        H1 = 79000;
        Gradient = 0;
    elseif h >= 90000 && h <= 100000
        T1 = 165.66;
        H1 = 90000;
        a = 0.004;
        Gradient = 1;
    elseif h > 100000
        disp('Please enter an altitude smaller than 100,000 meters')
    elseif h <= 0
        disp('Please enter a number larger than 0')
    else 
        disp('You broke my code :( (or dont enter a letter it wants a number)')
    end
    
    if Gradient == 1
        T = T1 + a*(h-H1);
        P = P1*(T/T1)^((-9.81)/(287.058*a));
        rho = rho1*(T/T1)^(-((9.81/287.058*a)+1));
    end
    
    if Gradient == 0
        T = T1 + a*(h-H1);
        P = P1*exp(-((9.81/(287.058*T)*(h-H1))));
        rho = rho1*exp(-((9.81/(287.058*T)*(h-H1))));
    end


end