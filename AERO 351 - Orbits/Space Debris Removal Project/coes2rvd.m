function [R1,V1] = coes2rvd(a,ecc,inc,RAAN,ArgP,theta,mu)
    %COES2RV Summary of this function goes here
    %   Detailed explanation goes here
    h = (mu*(a*(1-ecc^2)))^(1/2);
    
    R = [(((h^2)/mu)*(1/(1+ecc*cosd(theta))))*cosd(theta);...
        (((h^2)/mu)*(1/(1+ecc*cosd(theta))))*sind(theta);0];
    V = [(mu/h)*-sind(theta);(mu/h)*(ecc+cosd(theta));0];
    
    [~,Q] = ECI2PERI(ArgP,inc,RAAN);
    
    R1 = Q'*R;
    V1 = Q'*V;
end