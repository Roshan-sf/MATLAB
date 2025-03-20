%% Roshan Jaiswal-Ferri
%Section - 02
%Aero 331 HW 3: 3/14/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Problem 2

y = 0.1143;
Vz1 = [1,10,100,1000,10000]; %N

for i = 1:length(Vz1)
    Vz = Vz1(i);
    syms q1 q2 q3 q4 q5 q6 d
    
    eq1 = q1-q5-q2 == -1.1963*Vz;
    eq2 = q2+q4-q1 == 1.1963*Vz;
    eq3 = q3+q6-q4 == 0.1603*Vz;
    eq4 = q5-q3-q6 == -0.1603*Vz;
    
    eq5 = (10*q1)+(13.887*q2)-(2.519*q3)-(7.932*q4)-(7.932*q5) == 0;
    eq6 = (10*q1)+(6.368*q2)+(10*q3)-(31.5*q6) == 0;
    
    eq7 = Vz*(d) == 0.126*q1 + 0.054*q3 + 0.08*q4 + 0.08*q5 + 0.0268*q6;
    
    soln = solve([eq1,eq2,eq3,eq4,eq5,eq6,eq7],q1,q2,q3,q4,q5,q6,d);
    
    q_1(i) = double(soln.q1);
    q_2(i) = double(soln.q2);
    q_3(i) = double(soln.q3);
    q_4(i) = double(soln.q4);
    q_5(i) = double(soln.q5);
    q_6(i) = double(soln.q6);
    d_(i) = double(soln.d);
end

for i = 1:length(Vz1)
    disp(['Results for Vz = ', num2str(Vz1(i)), ':']);
    disp(['q1: ', num2str(q_1(i))]);
    disp(['q2: ', num2str(q_2(i))]);
    disp(['q3: ', num2str(q_3(i))]);
    disp(['q4: ', num2str(q_4(i))]);
    disp(['q5: ', num2str(q_5(i))]);
    disp(['q6: ', num2str(q_6(i))]);
    disp(['d: ', num2str(d_(i))]);
    disp(' ');
end