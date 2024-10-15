%Roshan Jaiswal-Ferri
%Section - 01
%Aero 320 HW 2 - Rotation Matricies: 10/8/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: Finding Rotation Matrix:

rV = [6783; 3391; 1953]; %Position Vector km
vV = [-3.5; 4.39; 4.44]; %Vel Vector km/s

%Converting to F'LVLH

Zlvlh = -(rV/norm(rV));
Ylvlh = -(cross(rV,vV)/norm(cross(rV,vV)));
Xlvlh = cross(Ylvlh,Zlvlh);

%Creating Matrix with new vectors

Clvlh_eci = [Xlvlh, Ylvlh, Zlvlh];
disp(Clvlh_eci)

