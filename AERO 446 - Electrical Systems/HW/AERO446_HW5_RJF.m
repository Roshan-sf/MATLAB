%% Roshan Jaiswal-Ferri
%Section - 01
%Aero 446 HW5: 5/9/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Problem 1

c = physconst('LightSpeed'); % m/s
n = 0.6;
A = 1;
f = 300*10^6; %hz
gamma = c/f;

Gain = (4*pi*n*A)/gamma^2;
GaindB = 10*log10(Gain);

disp(['Gain in dB: ', num2str(GaindB)])
disp(' ')

%% Problem 2

% 1) Radiation intensity of an antenna in a given direction

% 2) Operating frequency of antenna in hz, can be found as f=c/g

% 3) Used to increase signal level for better reception, should be designed
% to produce as little noise as possible and be close to source

% 4) Noise figure cannot be negative, NF = SNRin/SNRout, where
% SNR = signal power / noise power. SNRin is always larger than SNRout (if
% SNRout is larger theres a different larger problem at hand), because of
% this SNR >= 1 which means log(SNR) > 0.

% 5) GdB = log10(G) -> G = 10^(GdB/10) -> 10^0.1 = 1.0233 which is 2.33%

%% Problem 3

% Gains (in dB) and Noise Figures (in dB)
gain_dB = [30, 20, 13];       % A, B, C
gain = log2Lin(gain_dB);
NF_dB = [3, 2, 1.5];        % A, B, C
NF = log2Lin(NF_dB);

% Convert to linear scale

configG = [
    gain(1) gain(2) gain(3);   % ABC
    gain(1) gain(3) gain(2);   % ACB
    gain(2) gain(1) gain(3);   % BAC
    gain(2) gain(3) gain(1);   % BCA
    gain(3) gain(1) gain(2)    % CAB
];

configF = [
    NF(1) NF(2) NF(3);   % ABC
    NF(1) NF(3) NF(2);   % ACB
    NF(2) NF(1) NF(3);   % BAC
    NF(2) NF(3) NF(1);   % BCA
    NF(3) NF(1) NF(2)    % CAB
];

Gtotal = gain(1)*gain(2)*gain(3); %Same for all
GtotaldB = 10*log10(Gtotal);

%using Friis Formula:
for i = 1:5
    G1 = configG(i,1);
    G2 = configG(i,2);
    G3 = configG(i,3);

    F1 = configF(i,1);
    F2 = configF(i,2);
    F3 = configF(i,3);

    Ftotal(i) = F1 + (F2 - 1)/G1 + (F3 - 1)/(G1 * G2);
end

disp(['Total Gain—same for all (linear): ', num2str(Gtotal)])
disp(['Total Gain—same for all (dB): ', num2str(GtotaldB)])
disp(' ')
disp('Total NF using Friis Formula:')
disp(['Config: ABC: ', num2str(Ftotal(1))])
disp(['Config: ACB: ', num2str(Ftotal(2))])
disp(['Config: BAC: ', num2str(Ftotal(3))])
disp(['Config: BCA: ', num2str(Ftotal(4))])
disp(['Config: CAB: ', num2str(Ftotal(5))])

%% Problem 5

% A+ would only be worth it if the system is extremely sensitive to noise.
% a 20dB drop in performance is huge, way larger than a NF drop from 3 to
% 0.7dB. Using A+ total gain would drop from 1,995,000 to 19,950 
% (63dB to 42.9 dB), and NF would only drop from 1.9959 to 1.2612 (3dB to
% 1.01dB) when in the same configuration of ABC.


%% Functions

function [out] = log2Lin(in)
    out = 10.^(in/10); 
end