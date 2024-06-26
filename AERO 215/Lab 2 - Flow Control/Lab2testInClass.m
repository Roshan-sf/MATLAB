%Roshan Jaiswal-Ferri
%Aero 215 Lab 2: 09/26/23

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%variables: 
%vel = velocity

%Constants

airDensity = 1.225; %kg/m^3

%run the dynamic pressure equation from 0 to 250

for vel = 0:10:250; %10 in the middle means increments in 10

    dynPres = .5*airDensity*vel^2;

    disp(vel);
    disp(dynPres);
    
end

randNum = rand;

%% part2

if randNum < 0.33;

    disp('too small')

elseif randNum >= 0.33 && randNum < 0.666;

    disp('just right')
else
    disp('too big')
end




% for i = 1:10; run from one to ten  
%     disp(i);
% end;

% while i < 10;
%     disp(num2str(i));
%     i = i+1;
% end

% h1 = 3500; %meters
% 
% if h > 0 && h < 11000 %meters
% 
%     %execute code
% elseif h <= 11000 && h < 25000
% 
%     %alternate code, have as many elseif as you want
% 
% else 
% 
%     %else would be the final else statement. 
% 
% end