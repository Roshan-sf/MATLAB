%Roshan Jaiswal-Ferri
%Aero 215 Lab 3: 09/28/23

%%

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Part 1

%variables for part 1: 
g = 0;
x = 1;
y = 1;
i = 0;
j = 0;
k = 0;
row1Col = 0;
row2Col = 0;
row3Col = 0;



for row = 1:3
    for col = 1:4
        randNum = rand;
        array_loop(row,col) = randNum;       
    end
end

disp(array_loop)

a = min(array_loop,[],2); %1 searches each column, 2 searches whole row
%disp(a)

 i = a(1,1); %sets i to min value from row 1
 j = a(2,1); %sets j to min value from row 2
 k = a(3,1); %sets k to min value from row 3
 
 %disp([num2str(i), num2str(j), num2str(k)])

while row1Col == 0
   if array_loop(1,y) ~= i
    y = y + 1;
   
   else 
    row1Col = y; 
    %disp(row1Col)
   end
end

y = 1;

while row2Col == 0
   if array_loop(2,y) ~= j
    y = y + 1;
   
   else 
    row2Col = y; 
    %disp(row2Col)
   end
end

y = 1;

while row3Col == 0
   if array_loop(3,y) ~= k
    y = y + 1;
   
   else 
    row3Col = y; 
    %disp(row3Col)
   end
end

disp(['The minimum Number from row one: (1,',num2str(row1Col),') is ', num2str(i)]);
disp(['The minimum Number from row Two: (2,',num2str(row2Col),') is ', num2str(j)]);
disp(['The minimum Number from row Three: (3,',num2str(row3Col),') is ', num2str(k)]);


%% Part 2

%variables for part 2:
Prius.color = 'white';
Prius.year = 2013;
Prius.make = 'toyota';
Prius.body = 'hatchback';
Prius.power = 134;
Prius.weight = 3165;

Civic.color = 'white';
Civic.year = 2015;
Civic.make = 'honda';
Civic.body = 'sedan';
Civic.power = 143;
Civic.weight = 2800;


priRatio = Prius.power / Prius.weight;
civRatio = Civic.power / Civic.weight;

if priRatio < civRatio
    disp(['The Civic has a larger power to weight ratio of: ', num2str(civRatio)]);
elseif priRatio == civRatio
    disp(['The power to weight ratio of each car is equal. Civic:',num2str(civRatio),' Prius:',num2str(priRatio)])
else
    disp(['The Prius has a larger power to weight ratio of: ', num2str(priRatio)]);
end

%% Part 3

t = 0:0.02:0.6;
y = 3*sin(2*pi*10*t);

figure
plot(t,y) %The reason the graph does not look like a sign wave is because
%matlab has been told to only process points every .02 t on a graph that 
%only goes to .6t. This means that Matlab only processes 30 (.6 /.02)points 
% to plot on the graph which makes it look like straight lines when they 
%are connected instead of a curve. If you increase this ratio (either by
%increasing the .6 or decreasing the .02) it will look a lot more lik a
%wave.

%% Part 4

t = 0:0.001:100;
y = 3*sin(2*pi*10*t);

figure
plot(t,y)%The reason the graph does not look like a sign wave is the
%opposite of last time. Last time there was not enough points (30), this
%time there are way too many: 100,000 (100 /.001)points. All these points 
% make the graph look like a giant block. If you decrease the ratio by
% either decreasing 100 to around 1 or increasing .001 it will look more
% recognizable as a sign wave.

%% Part 5



xValues = linspace(0, 2*pi, 25); %Creates an evenly spaces vector starting at 0 going to 2pi,
%and the last number determines how many total values are calculated if
%left blank it will default to 100, replace 2*pi with 360 and sin with sind
%to work with degrees instead.

yValues = sin(xValues);

figure; %Figure creates a popup wondow, which we can put a graph on
plot(xValues, yValues,'r',LineWidth=4) %'* showes the plot points with no line. two sashes (--) is a dashed line
xlabel('Angles [rad] from 0 to 2pi');
ylabel('Sin(x)')
title('Graph of Sin(x)')
grid on;
grid minor; %Adds more smaller grid lines












