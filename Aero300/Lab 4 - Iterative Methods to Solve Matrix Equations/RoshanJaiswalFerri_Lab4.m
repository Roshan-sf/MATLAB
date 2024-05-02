%Roshan Jaiswal-Ferri
%Section - 03
%Aero 300 Lab 4 - Iterative Methods to Solve Matrix Equations: 4/26/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
%clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: 
tol = 1e-6;

a = [3,-1,2];
b = [1,4,2];
c = [-1,3,5];

%A = [a;b;c];
B = [1;1;1];
%d = diag(A);
y = randi([2 9]);
z = randi([2 9]);
A = randi([-10 10],y,z);
s = size(A);
isDom = 0; %Diagonal Dominance 0 is dominant 1 is not

%% PART 2: Testing for Diagonal Dominance

for i = 1:s(1,:) %For loop testing for diagonal dominance
    x=0;
    if j > s(1,2) || i > s(1,1) %stopping the loop if the matrix is mxn and not nxn and it is looping past an existing column or row
        break
    end
    for j = 1:s(:,1)
        x = x + A(i,j);
    end
    x = x - A(i,i);
    if A(i,i) < x %if any of the rows are not dominant set isDom to not dom
        isDom = 1; 
    end
end

if isDom == 1
    disp('Warning: May not Converge')
elseif isDom == 0
    disp('Matrix Will converge')
else
    disp('error (isDom)')
end

