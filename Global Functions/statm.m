function [T, P, rho] = statm(h)
%Function for finding T P and rho up to 10km
%   Input is in meters, ordered output:
%   [T, P, rho] = statm(h)
%   Temp in Kelvin
%   Pressure in KPa
%   Density in kg/m^3

layer = 0;
g = 9.81; % is g
r = 287.058; %gas constant
P1 = 101.325; %Pressure at sea level in kPa
rho1 = 1.2250; %density in kg/m^3
T1 = 288.16; %Initial temperature at sea level
H1 = 11000; %end of first gradient layer
H2 = 25000; %end of first isothermic layer
H3 = 47000; %end of second gradient layer
H4 = 53000; %end of second isothermic layer
H5 = 79000; %end of third gradient layer
H6 = 90000; %end of third isothermic layer
H7 = 100000; %end of fourth gradient layer
a1 = -.0065; %lapse rate for first gradient layer
a2 = .003; %lapse rate for second gradient layer
a3 = -.0045; %lapse rate for third gradient layer
a4 = .004; %lapse rate for fourth gradient layer

 if h <= 11000 && h >= 0
        layer = 1;
    elseif h > 11000 && h <= 25000
        layer = 2;
    elseif h > 25000 && h <= 47000
        layer = 3;
    elseif h > 47000 && h <= 53000
        layer = 4;
    elseif h > 53000 && h <= 79000
        layer = 5;
    elseif h > 79000 && h <= 90000
        layer = 6;
    elseif h > 90000 && h <= 100000
        layer = 7;
    elseif h > 100000
        eroor('Please enter an altitude smaller than 100,000 meters')
    elseif h <= 0
        error('Please enter a number larger than 0')
    else 
        error('You broke my code :( (or dont enter a letter it wants a number)')
 end

if layer >= 1
    if h >= 0 && h < 11000
        H1 = h;
    end
    a = a1;
    
    T = T1 + a*(H1-0);       
    P = P1*(T/T1)^((-g)/(r*a));  
    rho = rho1*(T/T1)^-((g/(r*a)+1));

    T1 = T;
    P1 = P;
    rho1 = rho;
end

if layer >= 2
    if h >= 11000 && h < 25000
        H2 = h;
    end
   
    T = T1;
    P = P1*exp(-((g/(r*T1)*(H2-H1))));
    rho = rho1*exp(-((g/(r*T)*(H2-H1))));

    T1 = T;
    P1 = P;
    rho1 = rho;
end

if layer >= 3
    if h >= 25000 && h < 47000
        H3 = h;
    end

    a = a2;
    T = T1 + a*(H3-H2);       
    P = P1*(T/T1)^((-g)/(r*a));  
    rho = rho1*(T/T1)^-((g/(r*a)+1));

    T1 = T;
    P1 = P;
    rho1 = rho;
end

if layer >= 4
    if h >= 47000 && h < 53000
        H4 = h;
    end
    
    T = T1;
    P = P1*exp(-((g/(r*T1)*(H4-H3))));
    rho = rho1*exp(-((g/(r*T)*(H4-H3))));

    T1 = T;
    P1 = P;
    rho1 = rho;
end

if layer >= 5
    if h >= 53000 && h < 79000
        H5 = h;
    end

    a = a3;

    T = T1 + a*(H5-H4);       
    P = P1*(T/T1)^((-g)/(r*a));  
    rho = rho1*(T/T1)^-((g/(r*a)+1));

    T1 = T;
    P1 = P;
    rho1 = rho;
end

if layer >= 6
    if h >= 79000 && h < 90000
        H6 = h;
    end
    
    T = T1;
    P = P1*exp(-((g/(r*T1)*(H6-H5))));
    rho = rho1*exp(-((g/(r*T)*(H6-H5))));

    T1 = T;
    P1 = P;
    rho1 = rho;
end

if layer >= 7
    if h >= 90000 && h < 100000
        H7 = h;
    end

    a = a4;

    T = T1 + a*(H7-H6);       
    P = P1*(T/T1)^((-g)/(r*a)); 
    rho = rho1*(T/T1)^-((g/(r*a)+1));

end
