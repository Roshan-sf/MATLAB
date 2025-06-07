%% Roshan Jaiswal-Ferri
%Section - 01
%Aero 446 Quiz 2: 5/15/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command 

%% Variables
Re = 6378; %km
Se = 1366; %w/m^2
mu = 398600;
Altgeo = 35687;
Rgeo = 42164; %km

f = 40e9; %Frequency in Hz
b = 10e6; %50 mhz
c = 3e8; %Speed of light in m/s
lambda = c / f; % Wavelength in meters
eta = 0.55;
OBO = 1; %dB
Ll = 1;
Ptx = 10*log10(60); %W to dB

tatmosphere = 100; %km
lookangle = 5; %deg
Cn = 3; %from problem


%EIRP = Cn - (Gr/Ts) + Ls - 228.6 - (10*log10(b));


%% Problem 1
Gtx = 1; %not stated
EIRP = Ptx * Gtx * OBO *Ll

Ls = (tatmosphere/sind(lookangle))*0.065 %db Total path loss

Ttx = 290;
Trx = ((10^(3/10))-1)*Ttx; %Ground system noise temp in K
Ts = 10*log10(Ttx + Trx)

syms Gr
eq = Cn == EIRP + (Gr/Ts) - Ls + 228.6 - (10*log10(b));
soln = solve(eq,Gr);
Gr = double(soln) %ground system gain for C/N of 3 dB

G = (10^(20/10));
syms A
eq2 = G == 4*pi*A*eta/lambda^2;
soln2 = solve(eq2,A);
A = double(soln2);

d = sqrt(2*pi*A)

%% Problem 2

% False
% False
% False
% False
% True

%% Problem 3

Codebook = [00 	00000;
    01 	00101;
    10 	10111;
    11 	01111;];

Data = Codebook(:,1);
Words = Codebook(:,2);

Code1 = '00000';
Code2 = '00101';
Code3 = '10111';
Code4 = '01111';



Data1 = '00101';
Data2 = '10101';
Data3 = '10111';
Data4 = '11111';
Data5 = '01111';
Data6 = '00001';

check2 = check(Data2,Code1,Code2,Code3,Code4);
check4 = check(Data4,Code1,Code2,Code3,Code4);
check6 = check(Data6,Code1,Code2,Code3,Code4);


% almost done all possible patterns are displayed they just need to be 
% sorted for the minimum number(s) and put together

[small, idx] = min(check2);
for i = 1:length(check2)
    if i ~= idx && check2(i) == small
        idx2 = i;
    end
end

[small2, idx3] = min(check4);
for i = 1:length(check4)
    if i ~= idx3 && check4(i) == small
        idx4 = i;
    end
end

[small3, idx5] = min(check6);
for i = 1:length(check6)
    if i ~= idx5 && check6(i) == small
        idx5 = i;
    end
end

% Found all small numbers and their position corresponding to which code,
% there is only a max of two double minimums for each check so I dont need
% more indexes











%% Functions

function out = check(data,code1,code2,code3,code4)
ham1 = 0;
ham2 = 0;
ham3 = 0;
ham4 = 0;
    for i = 1:length(data)
      digit_char = data(i); % Extract each character (digit) from the character array
      check_dig1 = code1(i);
      check_dig2 = code2(i);
      check_dig3 = code3(i);
      check_dig4 = code4(i);

      if strcmp(check_dig1,digit_char) == 0
          ham1 = ham1 + 1;
      end

      if strcmp(check_dig2,digit_char) == 0
          ham2 = ham2 + 1;
      end     

      if strcmp(check_dig3,digit_char) == 0
          ham3 = ham3 + 1;
      end

      if strcmp(check_dig4,digit_char) == 0
          ham4 = ham4 + 1;
      end

    end
      ham = [ham1,ham2,ham3,ham4];
      disp(ham)
      out = ham;

      %[~, out] = min(ham); %finds the position of smallest ham
end

