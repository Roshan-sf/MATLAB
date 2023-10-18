function [V, D, P] = hw2_dragPower_Jaiswal_FerriRoshan(v_max, aircraft)

%Variables
R = 287.058; %Gas constant
v_maxft = v_max*1.68781; %1.68781 is the conversion ratio for knots to ft/s
W = aircraft.WS*aircraft.S; %Calculating weight in lbf

%Breaking out structures to variabeles so equations are not so long
WS = aircraft.WS; %in lbf/ft^2
AR = aircraft.AR;
CD0 = aircraft.CD0;
e = aircraft.e;
CLM = aircraft.CLM;
S = aircraft.S; %ft^2
rho = aircraft.rho;



V_stall = sqrt((W)/((CLM)*(.5)*(rho)*(S))); %in ft/s

%disp(['vstall:', num2str(V_stall)]) was used for debugging purposes

V = linspace(V_stall,v_maxft,1000);

    for i = 1:length(V)
        D(i) = ((CD0)*(S)*(.5)*(rho)*(V(i)^2)+((W^2)/((.5)*(rho)*(V(i)^2)*(pi)*(AR)*(S)*(e)))); %in lbf
        P(i) = ((D(i))*(V(i)))/(550); %in hp
    end

    
V = V*0.592484;
    
end