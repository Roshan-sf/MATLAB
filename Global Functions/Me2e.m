function [Esoln] = Me2e(Me,ecc)
    %ME2E Summary of this function goes here
    %   Solves For Eccentric Anamoly (EA)
    %   Radians!
    
    syms E
    eq = Me == E - ecc*sin(E);
    Esoln = vpasolve(eq,E); %bounds ,[]
    Esoln = double(Esoln);

end