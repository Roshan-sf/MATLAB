function [r,v,DCM_ECI_PERI,DCM_PERI_ECI] = coes2rv(ecc,TA,inc,RAAN,w,mu,radius)
    %INPUT RADIANS!!!!
    %r1 = [cos()]

    wd = rad2deg(w);
    incd = rad2deg(inc);
    RAANd = rad2deg(RAAN);
    r = 1;
    v = 1;

    %h = sqrt(mu*(1+ecc)*radius);
    
    %r = (h^2/mu)*(1/(1+ecc*cos(TA)))

    DCM_ECI_PERI1 = rotz(wd)*rotx(incd)*rotz(RAANd);
    

    DCM_PERI_ECI = DCM_ECI_PERI';
end