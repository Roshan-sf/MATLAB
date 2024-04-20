%Roshan Jaiswal-Ferri
%Section - 03
%Aero 300 Lab 3 - Advance Data Loading and Plotting: 4/19/24

%% Clear Workspace

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: Using Bracketing Method


fx = @(x) x^3+1.0142*x^2-19.3629*x+15.8398;

d = 3; %Degree of polynomial
s = 1; %Step size
g = -6; %inital guess


%M = roots()

%% Bracketing Function

%rootMatrix = roots(initial_guess, step_size, expected_num_roots, function)

for row = 1:d %two for loops creating a 3x4 matrix
    for col = 1:2
        M(row,col) = 0;       
    end
end

lBracket = g;

% for i = 1:d
%     hBracket = lBracket+s;
%     h = fx(hBracket);
%     l = fx(lBracket);
%     if h*l<0 %if it is <0 then ans is negative and there is root in bracket
%         while done == 0 
%             hBracket2 = lBracket + e-6;
%             h = fx(hBracket2);
%             l = fx(lBracket);
%             if h*l<0
% 
% 
%         end
%     elseif h*l==0
%         %code
%     else
%         %code
%     end
    





%end







