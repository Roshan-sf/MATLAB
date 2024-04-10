%Section - 03
%Aero 300 Lec 3 - Bisection Method to find Roots: 4/10/24

%%

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: Bisection Loop

f = @(X)exp(x-2)+x.^3-x; %@(x) is a function handle
a = -2;
b = 2;

x = a:.01:b;
y = f(x);
plot(x,Y)

while ((b-a)/2) > e-5

end


%% Ex 2 In Class FPI

%% Newton-Raphson

h = @(x)x^3+x-1;
tic
[XVEC, XDIF, FX, NIT] = mynewton(0, 100, @(x)x^3+x-1, )

fzero() %used to find roots


format long
