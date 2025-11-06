%% Load Data:
data = readcell("495B434E.TXT");
% Units are in 1e-6
% Pre Opening:
hoop0 = data(7:19,4);
hoopStrain = cell2mat(hoop0); % 1e-6
avgHoopStrain0 = mean(hoopStrain);
hoopC0 = data(7:19,2);
hoopCStrain = cell2mat(hoopC0); % 1e-6
avgHoopCStrain0 = mean(hoopCStrain);
long0 = data(7:19,5);
longStrain0 = cell2mat(long0); % 1e-6
avgLongStrain0 = mean(longStrain0);
longC0 = data(7:19,3);
longCStrain0 = cell2mat(longC0); % 1e-6
avgLongCStrain0 = mean(longCStrain0);
% After opening:
hoopF = data(20:78,4);
hoopStrainF = cell2mat(hoopF); % 1e-6
avgHoopStrainF = mean(hoopStrainF);
longF = data(20:78,5);
longStrainF = cell2mat(longF); % 1e-6
avgLongStrainF = mean(longStrainF);
hoopCF = data(49:78,2);
hoopCStrainF = cell2mat(hoopCF); % 1e-6
avgHoopCStrainF = mean(hoopCStrainF);
longCF = data(49:78,3);
longCStrainF = cell2mat(longCF); % 1e-6
avgLongCStrainF = mean(longCStrainF);
% Can measurements:
jDiam = 65.85/1000; % m
tJ = [.118 .121 .130 .115 .134 .126 ]/1000; % m
avg_tJ = mean(tJ);
stdJ = std(tJ);
cDiam = 66.01/1000; % m
tC =[ .108 .155 .166 .130 .131 ]/1000;  % m
avg_tC = mean(tC);
stdC = std(tC);
% Pressure and Stress Calculations:
% Given Properties for alumminum can alloy 3004-H19
E = 69e9; % elastic modulus [Pa]
v = 0.35; % poisson's ratio
% strain measured:
% Juliet:
eLong_J  = (avgLongStrain0  - avgLongStrainF)  * 1e-6; % Pa
eHoop_J  = (avgHoopStrain0  - avgHoopStrainF)  * 1e-6; % Pa
% Cliodhna:
eLong_C  = (avgLongCStrain0 - avgLongCStrainF) * 1e-6; % Pa
eHoop_C  = (avgHoopCStrain0 - avgHoopCStrainF) * 1e-6; % Pa
% Finding internal pressure of can:
pJ_long = (4*avg_tJ*E*eLong_J) /(jDiam*(1 - 2*v)); % Pa;
pJ_hoop = (2*avg_tJ*E*eHoop_J) /(jDiam*(1 - v/2)); % Pa;
pC_long = (4*avg_tC*E*eLong_C) /(cDiam*(1 - 2*v)); % Pa;
pC_hoop = (2*avg_tC*E*eHoop_C) /(cDiam*(1 - v/2)); % Pa;
% Average internal presssure calcs with both hoop and long stresses
if (pJ_long < 0 && pJ_hoop < 0)
    pJ_long = -pJ_long;
    pJ_hoop = -pJ_hoop;
end
if (pC_long < 0 && pC_hoop < 0)
    pC_long = -pC_long;
    pC_hoop = -pC_hoop;
end
% internal pressure estimates in Pa
pJ = mean([pJ_long, pJ_hoop]) ;
pC = mean([pC_long, pC_hoop]) ;
pJ_long = pJ_long*1e-3; % KPa;
pJ_hoop = pJ_hoop*1e-3; % KPa;
pC_long = pC_long*1e-3; % KPa;
pC_hoop = pC_hoop*1e-3; % KPa;
% Finding hoop and long stresses:
sigmaHoopJ = (pJ * (jDiam/2)) / avg_tJ; % Pa
sigmaHoopC = (pC * (cDiam/2)) / avg_tC; % Pa
sigmaLongJ = (pJ * (jDiam/2)) / (2*avg_tJ); % Pa
sigmaLongC = (pC * (cDiam/2)) / (2*avg_tC); % Pa
stressRatioJ = sigmaHoopJ/sigmaLongJ; % Pa
stressRatioC = sigmaHoopC/sigmaLongC; % Pa
% Convert stresses from Pa to MPa
sigmaHoopJ = sigmaHoopJ * 1e-6;
sigmaHoopC = sigmaHoopC * 1e-6;
sigmaLongJ = sigmaLongJ * 1e-6;
sigmaLongC = sigmaLongC * 1e-6;
% Calculating Von Mises Stress:
% for thin walled pressure vessels "radial" stress = 0
sigma_r = 0;
% sigma1 = sigmaHoop;
% sigma2 = sigmaLong;
% thin walled simplified equation:
vonMisesJ = sqrt(sigmaHoopJ^2 + sigmaLongJ^2 - sigmaHoopJ*sigmaLongJ); % MPa
vonMisesC = sqrt(sigmaHoopC^2 + sigmaLongC^2 - sigmaHoopC*sigmaLongC); % MPa
% Calculating Factor of Safety:
yieldStress = 285; % MPa
fosJ = yieldStress/vonMisesJ;
fosC = yieldStress/vonMisesC;
disp("Juliet's Data:")
disp("average hoop Strain before opening can is " + avgHoopStrain0)
disp("average long Strain before opening can is " + avgLongStrain0)
disp("average hoop Strain after opening can is " + avgHoopStrainF)
disp("average long Strain after opening can is " + avgLongStrainF)
disp(" ")
disp("The diameter of Juliet's can is " + jDiam + " m")
disp("The average thickness of Juliet's can after opening is " + avg_tJ + " m")
disp("The standard deviation was " + stdJ)
disp(" ")
disp("The internal hoop pressure of Juliet's can is " + pJ_hoop + " KPa")
disp("The internal long pressure of Juliet's can is " + pJ_long + " KPa")
disp("The hoop stress of Juliet's can is " + sigmaHoopJ + " MPa")
disp("The long stress of Juliet's can is " + sigmaLongJ + " MPa")
disp("The stress ratio of Juliet's can is " + stressRatioJ)
disp(" ")
disp("Juliet's Von Mises Stress is " + vonMisesJ + " MPa")
disp("Juliet's Factor of Safety is " + fosJ)
disp(" ")
disp("Cliodhna's data:")
disp("average hoop Strain before opening can is " + avgHoopCStrain0)
disp("average long Strain before opening can is " + avgLongCStrain0)
disp("average hoop Strain after opening can is " + avgHoopCStrainF)
disp("average long Strain after opening can is " + avgLongCStrainF)
disp(" ")
disp("The diameter of Cliodhna's can is " + cDiam + " m")
disp("The average thickness of Cliodhna's can after opening is " + avg_tC + " m")
disp("The standard deviation was " + stdC)
disp(" ")
disp("The internal hoop pressure of Cliodhna's can is " + pC_hoop + " KPa")
disp("The internal long pressure of Cliodhna's can is " + pC_long + " KPa")
disp("The hoop stress of Cliodhna's can is " + sigmaHoopC + " MPa")
disp("The long stress of Cliodhna's can is " + sigmaLongC + " MPa")
disp("The stress ratio of Cliodhna's can is " + stressRatioC)
disp(" ")
disp("Cliodhna's Factor of Safety is " + fosC)
disp("Cliodhna's Von Mises Stress is " + vonMisesC + " MPa")
%disp("Our strain gauges were aligned perfectly when measured")