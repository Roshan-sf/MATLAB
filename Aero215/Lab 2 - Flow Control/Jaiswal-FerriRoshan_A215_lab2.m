%Roshan Jaiswal-Ferri
%Aero 215 Lab 2: 09/26/23

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%variables: 
%vel = velocity
%airDensity = is air density, currently set to desnity at sea level
%dynPres = dynamic pressure found by multiplying .5*airDensity*vel^2
%randNum = random number between 0 and 1


%Constants (Theres only one)

airDensity = 1.225; %kg/m^3

%run the dynamic pressure equation from 0 to 250

for vel = 0:10:250; %10 in the middle means increments in 10
                    %Calculating the dynamic air pressure using velocity
                    %in 10m/s^2 increments and displaying both

    dynPres = .5*airDensity*vel^2;

    disp(['Velocity is: ',num2str(vel),', Dynamic Pressure is:',num2str(dynPres)]);
    
end

%% part2

for i = 1:20; %compares a random # between 0-1, 
    % if its too small (0-.33), just right (.33-.66), too big (.66 up)
  
    randNum = rand;

    if randNum < 0.33;
        disp(['too small, randNum is: ',num2str(randNum),', Times Run: ', num2str(i)]);
    
    elseif randNum >= 0.33 && randNum < 0.666;

        disp(['Just Right, randNum is: ',num2str(randNum),', Times Run: ', num2str(i)]);
    else
        disp(['too big, randNum is: ',num2str(randNum),', Times Run: ', num2str(i)]);

    end
end


















