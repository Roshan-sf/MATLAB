clc
clear
close all
% All Groups Data
% Group 1
I_1 = [.01 .02 .03 .04 .05 .06 .07 .08 .09 .10 .11 0.12 .13 .14 .15 .16 .17 .18 .19 .2 .21];
V_1 = [.531 .497 .475 .454 .438 .415 .386 .374 .358 .334 .31 .287 .267 .242 .215 .189 .164 .136 .107 .09 .065];

Icon_1 = [.01 .02 .03 .04 .05 .06 .07 .08 .09];
Vcon_1 = [.482 .38 .324 .264 .218 .153 .108 .042 .044];

% Group 2
I_2 = [.01 .02 .03 .04 .05 .06 .07 .08 .09 .10 .11 0.12 .13 .14 .15 .16 .17];
V_2 = [.491 .455 .416 .373 .333 .282 .264 .209 .182 .166 .146 .123 .1 .081 .076 .054 .051];

Icon_2 = [.01 .02 .03 .04 .05 .06 .07 .08 .09 .10 .11 0.12 .13];
Vcon_2 = [.464 .426 .392 .353 .327 .291 .257 .218 .192 .152 .111 .062 .051];

% Group 3
I_3 = [.01 .02 .03 .04 .05 .06 .07 .08 .09 .10 .11 0.12 .13];
V_3 = [.471 .433 .423 .402 .389 .356 .336 .315 .301 .262 .226 .186 .14];

Icon_3 = [.01 .02 .03 .04 .05 .06 .07 .08 .09 .10 .11 0.12 .13];
Vcon_3 = [.484 .455 .425 .395 .372 .34 .307 .271 .228 .13 .094 .09 .057];

% Group 4

I_4 = [.01 .02 .03 .04 .05 .06 .07 .08 .09 .10 .11];
V_4 = [.472 .432 .391 .347 .312 .268 .222 .171 .117 .093 .053];
T_4 = [23.1 25.7 25.1 24.0 23.2 23.1 24.4 23.1 23.3 23.9 24.7]; 

Icon_4 = [.01 .02 .03 .04 .05 .06 .07 .08];
Vcon_4 = [.438 .378 .319 .261 .225 .163 .087 .043];
Tcon_4 = [23.5 23.0 24.0 24.1 25.2 24.4 23.7 23.2];

Iother_4 = [.01 .02 .03 .04 .05 .06 .07 .08 .09];
Vother_4 = [.479 .450 .420 .388 .363 .324 .278 .225 .170];

Iconother_4 = [.01 .02 .03 .04 .05 .06 .07 .08 .09 .10];
Vconother_4 = [.468 .433 .396 .360 .322 .281 .231 .199 .128 .095];


% Run Function for Each Group

% Group 1
[FF_1reg, eta_1reg, sigFF_1reg] = Lab2calcs(I_1,V_1)
[FF_1con, eta_1con, sigFF_1con] = Lab2calcs(Icon_1,Vcon_1)

% Group 2
[FF_2reg, eta_2reg, sigFF_2reg] = Lab2calcs(I_2,V_2)
[FF_2con, eta_2con, sigFF_2con] = Lab2calcs(Icon_2,Vcon_2)

% Group 3
[FF_3reg, eta_3reg, sigFF_3reg] = Lab2calcs(I_3,V_3)
[FF_3con, eta_3con, sigFF_3con] = Lab2calcs(Icon_3,Vcon_3)

% Group 4
[FF_4reg, eta_4reg, sigFF_4reg] = Lab2calcs(I_4,V_4)
[FF_4con, eta_4con, sigFF_4con] = Lab2calcs(Icon_4,Vcon_4)
[FF_4other, eta_4other, sigFF_4other] = Lab2calcs(Iother_4,Vother_4)
[FF_4conother, eta_4conother, sigFF_4conother] = Lab2calcs(Iconother_4,Vconother_4)


% Plot

figure(1)
plot(V_4,I_4,'-o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'LineWidth', 1.2, 'MarkerSize', 3)
xlabel("Voltage (V)")
ylabel("Current (A)")
title("Voltage vs Current for Clear Glass and Long View Factor")
hold on
errorbar(V_4,I_4,sigFF_4reg/100)
grid on

figure(2)
plot(Vcon_4,Icon_4,'-o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'LineWidth', 1.2, 'MarkerSize', 3)
xlabel("Voltage (V)")
ylabel("Current (A)")
title("Voltage vs Current for Contaminated Glass and Long View Factor")
hold on
errorbar(Vcon_4,Icon_4,sigFF_4con/100)
grid on

figure(3)
plot(Vconother_4,Iconother_4,'-o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'LineWidth', 1.2, 'MarkerSize', 3)
xlabel("Voltage (V)")
ylabel("Current (A)")
title("Voltage vs Current for Contaminated Glass")
hold on
errorbar(Vconother_4,Iconother_4,sigFF_4conother/100)
grid on

function [FillFactor, eta, sigFF] = Lab2calcs(I,V)

% Locate max values of V and I vectors
Vtot = max(V); % V
Itot = max(I); % I

% Calculate P tot
Ptot = Vtot*Itot; 

% Locate points on graph for Pmax
VtimesI = V.*I;
Pmax = max(VtimesI);

% Calculate Fill Factor
FillFactor = Pmax/Ptot;

DiameterLightBulb = 0.1778; % m
radiusLightBulb = DiameterLightBulb/2;

Plightbulb = 100; % W

A2 = 0.055*0.028;
A1 = pi*radiusLightBulb^2;
PIcirc = Plightbulb/(A1);
Pin = A2/A1 * PIcirc; % W/m^2

eta = (Pmax/(A2))/Pin * 100; % percent

% Error
[~, max_index] = max(VtimesI); 
n = max_index;
Imp = I(n);
Vmp = V(n);
Isc = I(end);
Voc = V(1);

sigImp = 0.01;
sigVmp = 0.001;
sigIsc = 0.01;
sigVoc = .001;

sigFF = sqrt((Vmp/(Isc*Voc)*sigImp)^2 + (Imp/(Isc*Voc)*sigVmp)^2 + (-(Imp*Vmp)/(Isc^2*Voc)*sigIsc)^2 + (-(Imp*Vmp)/(Isc*Voc^2)*sigVoc)^2);

end
