function [R1,V1,R2,V2] = coes2rvd(a,ecc,inc,RAAN,ArgP,theta,mu)
    %COES2RV The outputs are the same except transposed
    % for ease of use with 1x3 or 3x1 vectors 
    % (my old code used the first 2)
    %   Input COEs Get R & V
    
    h = (mu*(a*(1-ecc^2)))^(1/2);
    
    R = (h^2/mu)/(1+ecc*cosd(theta)) *[cosd(theta);sind(theta);0];
    V = (mu/h)*[-sind(theta);ecc+cosd(theta);0];
    
    [~,Q] = ECI2PERI(ArgP,inc,RAAN);
    
    R1 = Q*R;
    V1 = Q*V;
    
    R2 = R1';
    V2 = V1';
end
