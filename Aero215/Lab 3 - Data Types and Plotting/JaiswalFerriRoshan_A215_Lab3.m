%Roshan Jaiswal-Ferri
%Aero 215 Lab 3: 09/28/23

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%variables: 
g = 0;
x = 1;
y = 1;

for row = 1:3
    for col = 1:4
        randNum = rand;
        array_loop(row,col) = randNum;       
    end
end

disp(array_loop)

a = min(array_loop,[],2); %1 searches each column 2 searches whole row
disp(a)

i = a(1,1);
j = a(2,1);
k = a(3,1);

disp([num2str(i), num2str(j), num2str(k)])

while g = 0;
   if array_loop(x,:) == i;
    

   else


end

























