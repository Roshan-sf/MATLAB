function [V1,V2] = lambUVBi(R1,R2,dtime,Tm,mu,tol)
    %LAMBUVBI: Lambert Orbit Transfer With Universal Variable
    %   [V1,V2] = lambUVBi(R1,R2,dtime,Tm,mu,tol)
    %
    %   Given position 1 and 2, the time in between and mu
    %   Returns Velocity at positions 1 and 2
    %   Tm is 1 or -1, decides between long and short transfers
    %
    %   This function uses the bisection method to iterate

    R1n = norm(R1);
    R2n = norm(R2);
    Z = 0;
    C = 1/2;
    S = 1/6;
    Zu = 4*pi^2; %upper bound
    Zl = -4*pi^2; %lower bound
    dtl = 1; %change in time of loop (random guess)

    %deltaTheta = acos((dot(R1,R2)/(R1n*R2n))); %Alt eqs for A
    %A = sin(deltaTheta)*sqrt((R1n*R2n)/(1-cos(deltaTheta)));

    A = Tm*sqrt(R1n*R2n*(1+(dot(R1,R2)/(R1n*R2n))));
    
    while abs(dtl-dtime) > tol
        Y = R1n+R2n+(A*((Z*S-1)/sqrt(C)));        
        UV = sqrt(Y/C);
        dtl = (((UV^3)*S)/sqrt(mu))+((A*sqrt(Y))/sqrt(mu));

        if dtl < dtime
            Zl = Z; %reset zlower
        elseif dtl > dtime
            Zu = Z; %reset zupper
        end

        Z = 0.5*(Zu+Zl); %update z to midpoint
        [C,S] = stumpff(Z); %update stumpff c(z) s(z)

        f = 1-(((UV^2)/R1n)*C);
        g = dtl - (((UV^3)/sqrt(mu)) * S);
        fd = (sqrt(mu)/(R1n*R2n))*UV*((Z*S)-1);
        gd = 1-(((UV^2)/R2n)*C);
        %<3: f*gd - fd*g = 1

        %g = (1/mu)*(((Y/C)^(3/2))*S+(A*sqrt(Y)))-(1/mu)*((Y/C)^(3/2))*S;
        %g = A*sqrt(Y/mu);

        for i = 1:3
            V1(i) = (1/g)*(R2(i)-f*R1(i));
            V2(i) = (fd*R1(i))+(gd*V1(i));
        end
  
    end
end