%% Roshan Jaiswal-Ferri
%Aero 402 Homework 4: 10/27/25

%% Workspace Prep

warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Question 2: Non-Ideal Conditions

g2 = 1.4;
r2 = 297;
tc2 = 298;

g = 1.2;
R = 260; %J/kg-k
Tc = 1650; %combustion temp
tb = 10000;
F = 40; %N
g0 = 9.81; %at sea level

PcPe = linspace(1,100);
PePc = linspace(0,1);

CF = C_f(g,PePc, 1, 0); 
Cstar = 1009.945158;
mdotv = F./(Cstar.*CF);
Isp = F./(mdotv.*g0);

CF2 = C_f(g2,PePc, 1, 0); 
Cstar2 = 434; %m/s
mdotv2 = F./(Cstar2.*CF2);
Isp2 = F./(mdotv2.*g0);
eps2 = expansionRatio(g2,PcPe);

Vtank = linspace(0,1);
Ptank = 1./Vtank; %make all constants = 1 for simplicity

% non ideal

mdot = F/(Cstar*1.4); %1.4 is the C_f chosen from the knee of the plot
disp(num2str(mdot));
ISP = F/(mdot*g0);
disp(num2str(ISP));
tb = 40000/F;
mp = mdot*tb;
disp(num2str(mp));

% Plots

% % This first plot is not what the question asks for (that ones next) but
% % this makes much more sense, especially because as the Pe/Pc ratio
% % approaches zero the C_f approaches the ideal case
% figure
% grid on
% hold on
% plot(PePc,CF)
% title('C_F vs Pe/Pc')
% xlabel('C_F')
% ylabel('Pe/Pc') 

figure
grid on
hold on
plot(PcPe, C_f(g,1./PcPe, 1, 0))
plot(PcPe, C_f(g2,1./PcPe, 1, 0))
title('C_F vs Pc/Pe')
xlabel('C_F')
ylabel('Pc/Pe')
legend('Hydrazine','Cold Gas')

figure
grid on
hold on
plot(Isp,eps)
plot(Isp2,eps2)
title('Isp vs \epsilon')
xlabel('Isp (s)')
ylabel('\epsilon')
legend('Hydrazine','Cold Gas')

figure
grid on
hold on
plot(Vtank,Ptank)
title('Tank Volume vs Tank Pressure')
xlabel('Tank Volume (m^3)')
ylabel('Tank Pressure (Pa)')

figure
grid on
hold on
plot(PcPe,eps)
title('Ae/At (\epsilon) vs Pc')
xlabel('\epsilon')
ylabel('Pc (Pa)')

%% Question 3: 

mdot = 0.05; %nitrogen, h2o2, hydrazine

gamman2 = 1.4;
gammah2 = 1.25;
gammahy = 1.2;

Tcn2 = 300;
Tch2 = 1300;
Tchy = 1650;

Rn2 = 297;
Rh2 = 290;
Rhy = 260;

epsilon = 20;

[PcPen2, cstarn2, cfn2, Ispn2, thrustn2] = performance(gamman2, Rn2, Tcn2, epsilon, mdot);
[PcPeh2, cstarh2, cfh2, Isph2, thrusth2] = performance(gammah2, Rh2, Tch2, epsilon, mdot);
[PcPehy, cstarhy, cfhy, Isphy, thrusthy] = performance(gammahy, Rhy, Tchy, epsilon, mdot);

Propellant = {'Nitrogen (N2)'; 'Hydrogen (H2)'; 'Hydrazine (Hy)'};
Gamma = [gamman2; gammah2; gammahy];
Tc = [Tcn2; Tch2; Tchy];
R = [Rn2; Rh2; Rhy];
PcPe = [PcPen2; PcPeh2; PcPehy];
cstar = [cstarn2; cstarh2; cstarhy];
cf = [cfn2; cfh2; cfhy];
Isp = [Ispn2; Isph2; Isphy];
Thrust = [thrustn2; thrusth2; thrusthy];

results = table(Propellant, Gamma, Tc, R, PcPe, cstar, cf, Isp, Thrust);

disp('-----------------------------------------------');
disp(' Propellant Performance Comparison ');
disp('-----------------------------------------------');
disp(results);

%c = cStar(1.4,297,298)

%% Question 4:

 % a) CEA Printout:
 %  *******************************************************************************
 % 
 %         NASA-GLENN CHEMICAL EQUILIBRIUM PROGRAM CEA2, FEBRUARY 5, 2004
 %                   BY  BONNIE MCBRIDE AND SANFORD GORDON
 %      REFS: NASA RP-1311, PART I, 1994 AND NASA RP-1311, PART II, 1996
 % 
 % *******************************************************************************
 % 
 % 
 % 
 % 
 % ### CEA analysis performed on Tue 28-Oct-2025 01:12:05
 % 
 % # Problem Type: "Rocket" (Infinite Area Combustor)
 % 
 % prob case = _______________8869 ro equilibrium
 % 
 % # Pressure (1 value):
 % p,atm= 2
 % # Supersonic Area Ratio (1 value):
 % supar= 20
 % 
 % # You selected the following reactants:
 % reac
 % name H2O2(L)           wt%=100.0000
 % 
 % # You selected these options for output:
 % # short version of output
 % output short
 % # Proportions of any products will be expressed as Mass Fractions.
 % output massf
 % # Heat will be expressed as siunits
 % output siunits
 % 
 % # Input prepared by this script:/var/www/sites/cearun/cgi-bin/CEARUN/prepareInpu
 % tFile.cgi
 % 
 % ### IMPORTANT:  The following line is the end of your CEA input file!
 % end
 % 
 % 
 % 
 % 
 % 
 %              THEORETICAL ROCKET PERFORMANCE ASSUMING EQUILIBRIUM
 % 
 %           COMPOSITION DURING EXPANSION FROM INFINITE AREA COMBUSTOR
 % 
 % Pin =    29.4 PSIA
 % CASE = _______________
 % 
 %             REACTANT                    WT FRACTION      ENERGY      TEMP
 %                                          (SEE NOTE)     KJ/KG-MOL      K  
 % NAME        H2O2(L)                      1.0000000   -187780.000    272.740
 % 
 % O/F=    0.00000  %FUEL=  0.000000  R,EQ.RATIO= 0.500000  PHI,EQ.RATIO= 0.000000
 % 
 %                 CHAMBER   THROAT     EXIT
 % Pinf/P            1.0000   1.8079   284.05
 % P, BAR            2.0265   1.1209  0.00713
 % T, K             1274.54  1130.83   347.32
 % RHO, KG/CU M    4.3364-1 2.7034-1 5.6024-3
 % H, KJ/KG        -5520.56 -5781.47 -7039.39
 % U, KJ/KG        -5987.88 -6196.10 -7166.74
 % G, KJ/KG        -19324.1 -18028.6 -10800.9
 % S, KJ/(KG)(K)    10.8302  10.8302  10.8302
 % 
 % M, (1/n)          22.676   22.676   22.676
 % (dLV/dLP)t      -1.00000 -1.00000 -1.00000
 % (dLV/dLT)p        1.0001   1.0000   1.0000
 % Cp, KJ/(KG)(K)    1.8460   1.7848   1.4318
 % GAMMAs            1.2479   1.2586   1.3442
 % SON VEL,M/SEC      763.7    722.4    413.7
 % MACH NUMBER        0.000    1.000    4.213
 % 
 % PERFORMANCE PARAMETERS
 % 
 % Ae/At                      1.0000   20.000
 % CSTAR, M/SEC               1037.7   1037.7
 % CF                         0.6961   1.6796
 % Ivac, M/SEC                1296.4   1816.0
 % Isp, M/SEC                  722.4   1742.9
 % 
 % 
 % MASS FRACTIONS
 % 
 % H2O              0.52962  0.52963  0.52963
 % *OH              0.00002  0.00000  0.00000
 % *O2              0.47036  0.47037  0.47037
 % 
 %  * THERMODYNAMIC PROPERTIES FITTED TO 20000.K
 % 
 % NOTE. WEIGHT FRACTION OF FUEL IN TOTAL FUELS AND OF OXIDANT IN TOTAL OXIDANTS

% b)
% The results partially match my results from the previous question. It is
% possible they mismatch because CEA specifically calculates and accounts
% for a 2atm input pressure where my calcs do not. CEA also includes more
% advanced/specific chemistry that my method ignores. But the Cf values are
% extremely close (i got 1.78 and CEA got 1.67). The Cstar values are also
% close, i got 933 and cea got 1037.

% c)
 %  *******************************************************************************
 % 
 %         NASA-GLENN CHEMICAL EQUILIBRIUM PROGRAM CEA2, FEBRUARY 5, 2004
 %                   BY  BONNIE MCBRIDE AND SANFORD GORDON
 %      REFS: NASA RP-1311, PART I, 1994 AND NASA RP-1311, PART II, 1996
 % 
 % *******************************************************************************
 % 
 % 
 % 
 % 
 % ### CEA analysis performed on Tue 28-Oct-2025 01:10:32
 % 
 % # Problem Type: "Rocket" (Infinite Area Combustor)
 % 
 % prob case = _______________6230 ro equilibrium
 % 
 % # Pressure (1 value):
 % p,atm= 2
 % # Supersonic Area Ratio (1 value):
 % supar= 20
 % 
 % # You selected the following reactants:
 % reac
 % name CH6N2(L),MMH      wt%=100.0000
 % 
 % # You selected these options for output:
 % # short version of output
 % output short
 % # Proportions of any products will be expressed as Mass Fractions.
 % output massf
 % # Heat will be expressed as siunits
 % output siunits
 % 
 % # Input prepared by this script:/var/www/sites/cearun/cgi-bin/CEARUN/prepareInpu
 % tFile.cgi
 % 
 % ### IMPORTANT:  The following line is the end of your CEA input file!
 % end
 % 
 % 
 % 
 % 
 % 
 %              THEORETICAL ROCKET PERFORMANCE ASSUMING EQUILIBRIUM
 % 
 %           COMPOSITION DURING EXPANSION FROM INFINITE AREA COMBUSTOR
 % 
 % Pin =    29.4 PSIA
 % CASE = _______________
 % 
 %             REACTANT                    WT FRACTION      ENERGY      TEMP
 %                                          (SEE NOTE)     KJ/KG-MOL      K  
 % NAME        CH6N2(L),MMH                 1.0000000     54200.000    298.150
 % 
 % O/F=    0.00000  %FUEL=  0.000000  R,EQ.RATIO= 0.000000  PHI,EQ.RATIO= 0.000000
 % 
 %                 CHAMBER   THROAT     EXIT
 % Pinf/P            1.0000   1.7649   174.21
 % P, BAR            2.0265   1.1482  0.01163
 % T, K              961.67   893.65   548.83
 % RHO, KG/CU M    3.2447-1 2.0176-1 3.9155-3
 % H, KJ/KG         1176.43   837.74 -1071.71
 % U, KJ/KG          551.87   268.62 -1368.81
 % G, KJ/KG        -13604.2 -12897.5 -9507.06
 % S, KJ/(KG)(K)    15.3698  15.3698  15.3698
 % 
 % M, (1/n)          12.802   13.056   15.359
 % MW, MOL WT        10.973   11.351   15.359
 % (dLV/dLP)t      -1.06895 -1.07498 -1.00014
 % (dLV/dLT)p        1.7639   1.8870   1.0016
 % Cp, KJ/(KG)(K)    8.5440   9.6574   2.3681
 % GAMMAs            1.2013   1.1902   1.2973
 % SON VEL,M/SEC      866.2    823.0    620.8
 % MACH NUMBER        0.000    1.000    3.416
 % 
 % PERFORMANCE PARAMETERS
 % 
 % Ae/At                      1.0000   20.000
 % CSTAR, M/SEC               1220.4   1220.4
 % CF                         0.6744   1.7375
 % Ivac, M/SEC                1514.5   2260.6
 % Isp, M/SEC                  823.0   2120.4
 % 
 % 
 % MASS FRACTIONS
 % 
 % CH4              0.13924  0.16362  0.34821
 % C2H6             0.00001  0.00000  0.00000
 % *H2              0.09618  0.09006  0.04373
 % NH3              0.00052  0.00045  0.00016
 % *N2              0.60761  0.60766  0.60791
 % C(gr)            0.15644  0.13819  0.00000
 % 
 %  * THERMODYNAMIC PROPERTIES FITTED TO 20000.K
 % 
 % NOTE. WEIGHT FRACTION OF FUEL IN TOTAL FUELS AND OF OXIDANT IN TOTAL OXIDANTS

 % i)
 % MMH has a higher exhaust velocity and Isp than H2O2. This is because the
 % molecular weight of the products of MMH are much smaller than H2O2,
 % which then scales with exit velocity. Even though H2O2 burns hotter it
 % is too heavy to accelerate fast enough & exit vel matters more

 % ii)
 % hydrazine may perform better but its reactants are very toxic compared 
 % to H2O2. H2O2 can decompose into water vapor and oxygen and H2O2 can
 % also be cheaper. 



%% Functions

function CF = C_f(gamma, pe_pc, Ae_At, p0_pc)
CF = sqrt( (2.*gamma.^2./(gamma-1)) .* (2./(gamma+1)).^((gamma+1)./(gamma-1)) ...
           .* (1 - pe_pc.^((gamma-1)./gamma)) ) ...
     + (pe_pc - p0_pc) .* Ae_At;
end

function epsilon = expansionRatio(gamma, PcPe)

epsilon = ((2./(gamma+1)).^(1./(gamma-1))) .* (PcPe.^(1./gamma)) .* ...
          (((gamma+1)./(gamma-1)) .* (1 - PcPe.^((1-gamma)./gamma))).^(-0.5);

end

function c = cStar(gamma, R, Tc)
    top = sqrt(gamma*R*Tc);
    base = 2/(gamma + 1);
    bottomroot = sqrt( base^((gamma+1)/(gamma-1)) );
    
    c = top / (gamma * bottomroot);

end

function [PcPe, cstar, cf, Isp, thrust] = performance(gamma, R, Tc, epsilon, mdot)
    
    g = 9.81;
    syms PcPe
    eq1 = epsilon == ((2./(gamma+1)).^(1./(gamma-1))) .* (PcPe.^(1./gamma)) .* ...
          (((gamma+1)./(gamma-1)) .* (1 - PcPe.^((1-gamma)./gamma))).^(-0.5);
    sol = solve(eq1, PcPe);
    PcPe = double(sol);
    PcPe = max(real(PcPe(PcPe>1 & imag(PcPe)==0)));

    cstar = cStar(gamma, R, Tc);
    cf = C_f(gamma,1/PcPe, 20, 0 );
    thrust = mdot*cstar*cf;
    Isp = thrust/(mdot*g);

end
