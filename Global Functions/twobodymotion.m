function dstate = twobodymotion(time,state,mu) %dstate is derivitve of state
%FUNCTION put in descrip    

    %define vars
    x = state(1);
    y = state(2);
    z = state(3);
    dx = state(4); %vel
    dy = state(5); %vel
    dz = state(6); %vel
    
    %mag of pos vector
    r = norm([x y z]);
    
    %accel: !!eqs of motion!!
    ddx = -mu*x/r^3;
    ddy = -mu*y/r^3;
    ddz = -mu*z/r^3;
    
    dstate = [dx; dy; dz; ddx; ddy; ddz];

end