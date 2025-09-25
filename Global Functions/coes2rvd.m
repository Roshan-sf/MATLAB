function [R1,V1,RT,VT] = coes2rvd(a,ecc,inc,RAAN,ArgP,theta,mu)
    %[R1,V1,RT,VT] = coes2rvd(a,ecc,inc,RAAN,ArgP,theta,mu)
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
    
    RT = R1';
    VT = V1';
end
