%Roshan Jaiswal-Ferri
%Section - 
%Aero 321 Lab 2: 10/11/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: Data interp

dist = readmatrix('DistanceData1');
fastDist = readmatrix('DistanceDataFast.csv');
fastDist = fastDist(532:end,:);

expected = 10; %cm

[avg1, stand_dev1, Err1] = interp(dist, expected);
[avg2, stand_dev2, Err2] = interp(fastDist, expected);


[h, p, stats] = chi2gof(dist(:,2), 'CDF', {@normcdf, avg1, stand_dev1});
[h2, p2, stats2] = chi2gof(fastDist(:,2), 'CDF', {@normcdf, avg2, stand_dev2});


function [avg, stand_dev, Err] = interp(dataVec, expected)

    avg = mean(dataVec(:,2));
    stand_dev = std(dataVec(:,2));
    Err = abs(((avg-expected)/expected)*100);

end
