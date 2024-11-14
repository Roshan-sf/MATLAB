%Roshan Jaiswal-Ferri
%Section - 01
%Aero 320 HW 4 - Problem 2: 10/21/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1:

RPM300A0 = load('S2G1RPM300AOAN0.mat');
RPM300A2 = load('S2G2RPM300AOAN2.mat');
RPM300A6 = load('S2G3RPM300AOAN6.mat');
RPM300A8 = load('S3G1RPM300AOAN08.mat');
RPM300A10 = load('S3G2RPM300AOAN10.mat');
RPM300A12 = load('S3G3RPM300AOAN12.mat');
RPM200A0 = load('S4G1RPM200AOAN00.mat');
RPM200A2 = load('S4G2RPM200AOAN02.mat');
RPM200A6 = load('S4G3RPM200AOAN06.mat');
RPM200A8 = load('S5G1RPM200AOAN8.mat');
RPM200A10 = load('S5G2RPM200AOAN10.mat');
RPM200A12 = load('S5G3RPM200AOAN12.mat');
RPM150A0 = load('S6G1RPM150AOAN00.mat');
RPM150A2 = load('S6G2RPM150AOAN02.mat');
RPM150A6 = load('S6G3RPM150AOAN06.mat');
RPM150A8 = load('S7G1RPM150AOAN08.mat');
RPM150A10 = load('S7G2RPM150AOAN10.mat');
RPM150A12 = load('S7G3RPM150AOAN12.mat');

RPM300A0P = mean(RPM300A0.P);
RPM300A2P = mean(RPM300A2.P);
RPM300A6P = mean(RPM300A6.P);
RPM300A8P = mean(RPM300A8.P);
RPM300A10P = mean(RPM300A10.P);
RPM300A12P = mean(RPM300A12.P);
RPM200A0P = mean(RPM200A0.P);
RPM200A2P = mean(RPM200A2.P);
RPM200A6P = mean(RPM200A6.P);
RPM200A8P = mean(RPM200A8.P);
RPM200A10P = mean(RPM200A10.P);
RPM200A12P = mean(RPM200A12.P);
RPM150A0P = mean(RPM150A0.P);
RPM150A2P = mean(RPM150A2.P);
RPM150A6P = mean(RPM150A6.P);
RPM150A8P = mean(RPM150A8.P);
RPM150A10P = mean(RPM150A10.P);
RPM150A12P = mean(RPM150A12.P);

RPM300A0F = mean(RPM300A0.F);
RPM300A2F = mean(RPM300A2.F);
RPM300A6F = mean(RPM300A6.F);
RPM300A8F = mean(RPM300A8.F);
RPM300A10F = mean(RPM300A10.F);
RPM300A12F = mean(RPM300A12.F);
RPM200A0F = mean(RPM200A0.F);
RPM200A2F = mean(RPM200A2.F);
RPM200A6F = mean(RPM200A6.F);
RPM200A8F = mean(RPM200A8.F);
RPM200A10F = mean(RPM200A10.F);
RPM200A12F = mean(RPM200A12.F);
RPM150A0F = mean(RPM150A0.F);
RPM150A2F = mean(RPM150A2.F);
RPM150A6F = mean(RPM150A6.F);
RPM150A8F = mean(RPM150A8.F);
RPM150A10F = mean(RPM150A10.F);
RPM150A12F = mean(RPM150A12.F);


RPM300A0q = RPM300A0P(1)-RPM300A0P(2);
RPM300A2q = RPM300A2P(1)-RPM300A2P(2);
RPM300A6q = RPM300A6P(1)-RPM300A6P(2);
RPM300A8q = RPM300A8P(1)-RPM300A8P(2);
RPM300A10q = RPM300A10P(1)-RPM300A10P(2);
RPM300A12q = RPM300A12P(1)-RPM300A12P(2);
RPM200A0q = RPM200A0P(1)-RPM200A0P(2);
RPM200A2q = RPM200A2P(1)-RPM200A2P(2);
RPM200A6q = RPM200A6P(1) - RPM200A6P(2);
RPM200A8q = RPM200A8P(1) - RPM200A8P(2);
RPM200A10q = RPM200A10P(1) - RPM200A10P(2);
RPM200A12q = RPM200A12P(1) - RPM200A12P(2);
RPM150A0q = RPM150A0P(1) - RPM150A0P(2);
RPM150A2q = RPM150A2P(1) - RPM150A2P(2);
RPM150A6q = RPM150A6P(1) - RPM150A6P(2);
RPM150A8q = RPM150A8P(1) - RPM150A8P(2);
RPM150A10q = RPM150A10P(1) - RPM150A10P(2);
RPM150A12q = RPM150A12P(1) - RPM150A12P(2);

tare = load('TARE_20mps_neg30_pos20.mat');
AOA = linspace(-30,20,3673);
tareP = mean(tare.P);
q20ms = tareP(1)-tareP(2);

% figure(1)
% plot(AOA,tare.F(:,1))
% hold on
% plot(AOA,tare.F(:,2))
% hold on
% plot(AOA,tare.F(:,3))
% legend('Fx','Fy','Fz','location','best')

tare0x = -1.22;
tare0z = 0.088;
tare2x = -1.20;
tare2z = 0.072;
tare6x = -1.18;
tare6z = 0.054;
tare8x = -1.16;
tare8z = 0.058;
tare10x = -1.158;
tare10z = 0.0175;
tare12x = -1.144;
tare12z = 0.024;

tarex = [tare0x tare2x tare6x tare8x tare10x tare12x];
tarez = [tare0z tare2z tare6z tare8z tare10z tare12z];

ratiox = tarex./q20ms;
ratioz = tarez./q20ms;

rpm300q = mean([RPM300A0q RPM300A2q RPM300A6q RPM300A8q RPM300A10q RPM300A12q]);
rpm200q = mean([RPM200A0q RPM200A2q RPM200A6q RPM200A8q RPM200A10q RPM200A12q]);
rpm150q = mean([RPM150A0q RPM150A2q RPM150A6q RPM150A8q RPM150A10q RPM150A12q]);

Fsting300x = ratiox*rpm300q;
Fsting200x = ratiox*rpm200q;
Fsting150x = ratiox*rpm150q;
Fsting300z = ratioz*rpm300q;
Fsting200z = ratioz*rpm200q;
Fsting150z = ratioz*rpm150q;

forcex300 = [RPM300A0F(1) RPM300A2F(1) RPM300A6F(1) RPM300A8F(1) RPM300A10F(1) RPM300A12F(1)];
forcex200 = [RPM200A0F(1) RPM200A2F(1) RPM200A6F(1) RPM200A8F(1) RPM200A10F(1) RPM200A12F(1)];
forcex150 = [RPM150A0F(1) RPM150A2F(1) RPM150A6F(1) RPM150A8F(1) RPM150A10F(1) RPM150A12F(1)];
forcez300 = [RPM300A0F(3) RPM300A2F(3) RPM300A6F(3) RPM300A8F(3) RPM300A10F(3) RPM300A12F(3)];
forcez200 = [RPM200A0F(3) RPM200A2F(3) RPM200A6F(3) RPM200A8F(3) RPM200A10F(3) RPM200A12F(3)];
forcez150 = [RPM150A0F(3) RPM150A2F(3) RPM150A6F(3) RPM150A8F(3) RPM150A10F(3) RPM150A12F(3)];


D300 = 4.448*(forcex300 - Fsting300x);
D200 = 4.448*(forcex200 - Fsting200x);
D150 = 4.448*(forcex150 - Fsting150x);
L300 = 4.448*(forcez300 - Fsting300z);
L200 = 4.448*(forcez200 - Fsting200z);
L150 = 4.448*(forcez150 - Fsting150z);

Sref = 0.670*0.115;

CD300 = D300./(Sref*rpm300q);
CD200 = D200./(Sref*rpm200q);
CD150 = D150./(Sref*rpm150q);
CL300 = L300./(Sref*rpm300q);
CL200 = L200./(Sref*rpm200q);
CL150 = L150./(Sref*rpm150q);

r = 0.115;
rhoinf = 1.225;

uinf300 = sqrt(2*rpm300q/rhoinf);
uinf200 = sqrt(2*rpm200q/rhoinf);
uinf150 = sqrt(2*rpm150q/rhoinf);

muinf = 1.789e-5;

Re300 = rhoinf*uinf300*2*r/muinf;
Re200 = rhoinf*uinf200*2*r/muinf;
Re150 = rhoinf*uinf150*2*r/muinf;

AOA = [0 2 6 8 10 12];
AOAlinspace = linspace(0,12,100);

coefL150 = polyfit(AOA,CL150,3);
yl150 = polyval(coefL150,AOAlinspace);

coefL200 = polyfit(AOA,CL200,3);
yl200 = polyval(coefL200,AOAlinspace);

coefL300 = polyfit(AOA,CL300,3);
yl300 = polyval(coefL300,AOAlinspace);

coefD150 = polyfit(AOA,CD150,3);
yD150 = polyval(coefD150,AOAlinspace);

coefD200 = polyfit(AOA,CD200,3);
yD200 = polyval(coefD200,AOAlinspace);

coefD300 = polyfit(AOA,CD300,3);
yD300 = polyval(coefD300,AOAlinspace);

for i = 1:length(CL300)
CDL300(i) = CL300(i)/CD300(i);
CDL200(i) = CL200(i)/CD200(i);
CDL150(i) = CL150(i)/CD150(i);
end

%% Error Calc

sigCL150 = sqrt(((1/rpm150q*Sref)^2)+((L150./(Sref*(rpm150q^2))).^2));
sigCL200 = sqrt(((1/rpm200q*Sref)^2)+((L200./(Sref*(rpm200q^2))).^2));
sigCL300 = sqrt(((1/rpm300q*Sref)^2)+((L300./(Sref*(rpm300q^2))).^2));

sigCD150 = sqrt(((1/rpm150q*Sref)^2)+((D150./(Sref*(rpm150q^2))).^2));
sigCD200 = sqrt(((1/rpm200q*Sref)^2)+((D200./(Sref*(rpm200q^2))).^2));
sigCD300 = sqrt(((1/rpm300q*Sref)^2)+((D300./(Sref*(rpm300q^2))).^2));

sigClCd150 = sqrt((((1./CD150).*sigCL150).^2)+(((CL150./(CD150.^2)).*sigCD150).^2));
sigClCd200 = sqrt((((1./CD200).*sigCL200).^2)+(((CL200./(CD200.^2)).*sigCD200).^2));
sigClCd300 = sqrt((((1./CD300).*sigCL300).^2)+(((CL300./(CD300.^2)).*sigCD300).^2));

%% Figures

colors = lines(3);

figure(1)
hold on
p1 = plot(AOA, CD150, 'Color', colors(1,:));
errorbar(AOA, CD150, sigCL150, 'o', 'Color', colors(1,:), 'HandleVisibility', 'off')
p2 = plot(AOA, CD200, 'Color', colors(2,:));
errorbar(AOA, CD200, sigCL200, 'o', 'Color', colors(2,:), 'HandleVisibility', 'off')
p3 = plot(AOA, CD300, 'Color', colors(3,:));
errorbar(AOA, CD300, sigCL300, 'o', 'Color', colors(3,:), 'HandleVisibility', 'off')
xlabel('Angle of Attack (deg)')
ylabel('Coefficient of Drag')
title('AOA vs Coefficient of Drag NACA 4412')
legend([p1, p2, p3], 'Re = 56,140', 'Re = 81,890', 'Re = 130,010', 'location', 'northwest')

figure(2)
hold on
p1 = plot(AOA, CL150, 'Color', colors(1,:));
errorbar(AOA, CL150, sigCD150, 'o', 'Color', colors(1,:), 'HandleVisibility', 'off')
p2 = plot(AOA, CL200, 'Color', colors(2,:));
errorbar(AOA, CL200, sigCD200, 'o', 'Color', colors(2,:), 'HandleVisibility', 'off')
p3 = plot(AOA, CL300, 'Color', colors(3,:));
errorbar(AOA, CL300, sigCD300, 'o', 'Color', colors(3,:), 'HandleVisibility', 'off')
xlabel('Angle of Attack (deg)')
ylabel('Coefficient of Lift')
title('AOA vs Coefficient of Lift NACA 4412')
legend([p1, p2, p3], 'Re = 56,140', 'Re = 81,890', 'Re = 130,010', 'location', 'best')

figure(3)
hold on
p1 = plot(AOA, CDL150, 'Color', colors(1,:));
errorbar(AOA, CDL150, sigClCd150, 'o', 'Color', colors(1,:), 'HandleVisibility', 'off')
p2 = plot(AOA, CDL200, 'Color', colors(2,:));
errorbar(AOA, CDL200, sigClCd200, 'o', 'Color', colors(2,:), 'HandleVisibility', 'off')
p3 = plot(AOA, CDL300, 'Color', colors(3,:));
errorbar(AOA, CDL300, sigClCd300, 'o', 'Color', colors(3,:), 'HandleVisibility', 'off')
xlabel('Angle of Attack (deg)')
ylabel('Cl/Cd')
title('AOA vs Cl/Cd NACA 4412')
legend([p1, p2, p3], 'Re = 56,140', 'Re = 81,890', 'Re = 130,010', 'location', 'best')




