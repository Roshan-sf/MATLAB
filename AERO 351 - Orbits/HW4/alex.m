%% Problem 8.16
% Givens
Date1_816 = [2005, 8, 15, 0]; % Departure date
Date2_816 = [2006, 3, 15, 0]; % Martian arrival date
zE_816 = 190; % [km] Parking orbit Earth altitude
incE_816 = 52; % [degrees] Inclination of the orbit
zpM_816 = 300; % [km] Ending Martian altitude of periapse
TM_816 = 35*60*60; % [s] Period of the ending orbit at Mars
rE_816 = zE_816 + rEarth; % [km]
rpM_816 = zpM_816 + rMars; % [km]
% Get Julian date numbers
[~,~,Tstart_816] = julianDate(Date1_816(1),Date1_816(2),Date1_816(3),Date1_816(4));
[~,~,Tend_816] = julianDate(Date2_816(1),Date2_816(2),Date2_816(3),Date2_816(4));
TLamberts_816 = (Tend_816 - Tstart_816)*24*60*60; % [s] Time for Lamberts transfer
Tstart_816 = Tstart_816/36525; % Convert to century year
Tend_816 = Tend_816/36525; % Convert to century year
% Get planetary COES
[Earth_coes_816] = AERO351planetary_elements2(3,Tstart_816);
[Mars_coes_816] = AERO351planetary_elements2(4,Tend_816);
% Breakout useful elements
eccEarth_816 = Earth_coes_816(2); % [~] Eccentricity of Earth
incEarth_816 = Earth_coes_816(3); % [degrees] Inclination of Earth
RAANEarth_816 = Earth_coes_816(4); % [degrees] Right ascension of Earth
eccMars_816 = Mars_coes_816(2); % [~] Eccentricity of Mars
incMars_816 = Mars_coes_816(3); % [degrees] Inclination of Mars
RAANMars_816 = Mars_coes_816(4); % [degrees] Right ascension of Mars
% Calculate more orbital elements for the Earth and Mars
[omegaEarth_816,thetaEarth_816,rpEarth_816] = planetCOEs(Earth_coes_816);
[omegaMars_816,thetaMars_816,rpMars_816] = planetCOEs(Mars_coes_816);
% Get R and V vectors for Earth and Mars
[REarth_816,VEarth_816,~,~,~] = COEstoRV(muSun,rpEarth_816,eccEarth_816,incEarth_816,RAANEarth_816,omegaEarth_816,thetaEarth_816);
[RMars_816,VMars_816,~,~,~] = COEstoRV(muSun,rpMars_816,eccMars_816,incMars_816,RAANMars_816,omegaMars_816,thetaMars_816);
% Call the Lamberts function to run Lamberts calculation
[V1_816,V2_816] = Lamberts(REarth_816,RMars_816,TLamberts_816,muSun);
vInfEarth_816 = norm(V1_816 - VEarth_816); % [km/s] Required hyperbolic excess speed from the Earth
vInfMars_816 = norm(V2_816 - VMars_816); % [km/s] Arrival hyperbolic excess speed at Mars
% Calculate delta V from the Earth
vPark_816 = sqrt(muEarth/rE_816); % [km/s]
vpHypEarth_816 = sqrt(vInfEarth_816^2 + 2*(muEarth/rE_816)); % [km/s]
deltaV1_816 = vpHypEarth_816 - vPark_816; % [km/s] Departure delta V
% Calculate delta V to arrive at Mars
aEnd_816 = (TM_816*sqrt(muMars)/(2*pi))^(2/3); % [km] Semimajor axis of the ending orbit
eccEnd_816 = 1 - rpM_816/aEnd_816; % [~] Eccentricity of the ending orbit
hEnd_816 = sqrt(muMars*aEnd_816*(1-eccEnd_816^2)); % Angular momentum of the ending orbit
vpEnd_816 = hEnd_816/rpM_816; % [km/s] Velocity at periapse of the ending orbit
vpHypMars_816 = sqrt(vInfMars_816^2 + 2*(muMars/rpM_816)); % [km/s] Velocity of the hyperbolic orbit at periapse of the ending orbit
deltaV2_816 = vpHypMars_816 - vpEnd_816; % [km/s] Arrival delta V
deltaVTot_816 = deltaV1_816 + deltaV2_816; % [km/s] Total mission delta V
% Display results
disp('<<<<< Problem 8.16 >>>>>');
disp(' ');
disp(['Delta V: ',num2str(deltaVTot_816),' km/s']);
disp(' ');
% Heart checks:
% Both the Earth and Mars COEs from Dr. A's function
% are reasonable, R and V vectors for both
% Mars and Earth match internet results and are reasonable, the required
% hyperbolic excess speed from earth and the hyperbolic excess speed upon
% Martian arrival are both reasonable, the required delta V to get into the
% hyperbolic orbit from the LEO parking orbit around Earth is reasonable,
% the eccentricty for the ending orbit at Mars is reasonbale given the
% radius of periapse and the given period, the velocity at periapse of that
% orbit is reasonable given how eccentric it is, the delta V required
% to go from the hyperbolic arrival to that orbit is also reasonable.
has context menu