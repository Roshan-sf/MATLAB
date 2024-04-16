%Roshan Jaiswal-Ferri
%Section - 03
%Aero 300 Lec 4 - LU: 4/15/24

%%

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Class Notes: LU

H = [3,2,1];
N = [9,8,7];

M = [H;N]; %creates a matrix ; means go down a row , is next column

[L,U] = lu(M);

disp(M)

disp(L)
disp(U)