function [T, P, rho] = IsenCalc(mode, totalMolW, M, P0, T0, rho0)
    %units:
    %total molecular weight in g/mol
    %Pressure in Pascals
    %Temperature in kelvin
    %Density in kg/m^3

    Ru = 8.314; %J/mol-K
    Rspecific = Ru/(totalMolW*1e-3); %convert from g to kg
    
    if mode == 1
        f = 3/2;
    elseif mode == 2
        f = 5/2;
    elseif mode == 3
        f = 6/2;
    else
        error('Mode Error')
    end

    Cv = f*Rspecific;
    Cp = Cv+Rspecific;
    Gamma = Cp/Cv;

    T = ((1+((Gamma-1)/2)*M^2)/T0)^-1;
    P = (((T0/T)^(Gamma/(Gamma-1)))/P0)^-1;
    rho = (((P0/P)^(1/Gamma))/rho0)^-1;

end