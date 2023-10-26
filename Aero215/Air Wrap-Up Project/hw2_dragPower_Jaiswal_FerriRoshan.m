function [V, D, P] = hw2_dragPower_Jaiswal_FerriRoshan(input2, input, rho1)

%Variables
R = 287.058; %Gas constant
%W = input.WS*input.S; %Calculating weight in lbf

%Breaking out structures to variabeles so equations are not so long
AR = input.AR;
CD0 = input.CD0;
e = input.e;
CLM = input.CLM;
S = input.S; %ft^2
rho = rho1*0.00194032033; %Converting rho from kg/m^3 to slug/ft^3
W = input.EW + input2.PW*input.PS + input.MF;
v_maxft = input.VM*1.68781; %1.68781 is the conversion ratio for knots to ft/s

V_stall = sqrt((W)/((CLM)*(.5)*(rho)*(S))); %in ft/s

disp([num2str(W), num2str(rho1)])
%disp(['vstall:', num2str(V_stall)]) was used for debugging purposes

V = linspace(V_stall,v_maxft,1000);

    for i = 1:length(V)
        D(i) = ((CD0)*(S)*(.5)*(rho)*(V(i)^2)+((W^2)/((.5)*(rho)*(V(i)^2)*(pi)*(AR)*(S)*(e)))); %in lbf
        P(i) = ((D(i))*(V(i)))/(550); %in hp
    end

    
V = V*0.592484; %ft/s to kts
    
end