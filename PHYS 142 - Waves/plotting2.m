close all
clear all
clc

pressure = [100.11
    101.82
    101.92
    100.02
    100.06]; %kPa

volume = [216.583125
    211.6081875
    249.749375
    254.7243125
    213.2665]; %cm^3

figure()
plot(pressure,volume)
grid on
xlabel('Pressure (kPa)')
ylabel('Volume (cm^3)')
title('PV')


pressure = [100.11
    101.82
    101.92
    100.02
    100.11]; %kPa

volume = [216.583125
    211.6081875
    249.749375
    254.7243125
    216.583125]; %cm^3

area = 0.5 * abs( sum(volume(1:end-1).*pressure(2:end) - volume(2:end).*pressure(1:end-1)) );

fprintf('Enclosed Area: %.4f kPa·cm^3\n', area);

area_Joules = area * 0.001;
fprintf('Enclosed Area: %.6f J\n', area_Joules);
