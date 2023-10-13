function [V, D, P] = hw2_dragPower_Jaiswal_FerriRoshan(v_max, aircraft)

%Variables
rho = .0023769; %rhoSealevel in slugs/ft^3
R = 287.058; %Gas constant
v_maxft = v_max*1.68781; %1.68781 is the conversion ratio for knots to ft/s
weight = aircraft.WS*aircraft.S; %Calculating weight in lbf

D = ((aircraft.CD0)*(aircraft.S)*(.5)*(rho)*(v_maxft))+((weight^2)/((.5)*(rho)*(v_maxft^2)*(pi)*(aircraft.AR)*(aircraft.S)*(aircraft.e)));








end