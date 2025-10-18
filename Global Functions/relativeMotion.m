function dstate = relativeMotion(time,state,mu) 
%Use Column vectors!
%INPUTS: first 6 rows: [x y z dx dy dz...] in ECI, target properties
%CONTD: second 6 rows: [...x y z dx dy dz] in LVLH, relative to target
%OUTPUT: follows same convention

    %unpack for clarity (t for target c for chaser):
    tx0 = state(1); %pos
    ty0 = state(2);
    tz0 = state(3);
    tdx0 = state(4); %vel
    tdy0 = state(5);
    tdz0 = state(6);

    cx0 = state(7); %pos
    cy0 = state(8);
    cz0 = state(9);
    cdx0 = state(10); %vel
    cdy0 = state(11);
    cdz0 = state(12);
    
    rvect = [tx0 ty0 tz0]; %r and v vectors for chaser from chief
    vvect = [tdx0 tdy0 tdz0];

    rt = norm([tx0 ty0 tz0]); %r vector magnitudes
    rc = rt; %norm([cx0 cy0 cz0]);
    hc = norm(cross(rvect,vvect));

    %target
    tddx = -mu*tx0/rt^3;
    tddy = -mu*ty0/rt^3;
    tddz = -mu*tz0/rt^3;

    dstate_t = [tdx0; tdy0; tdz0; tddx; tddy; tddz];

    %chaser
    cddx = ((2*mu/rc^3)+(hc^2/rc^4))*cx0 - 2*(dot(vvect,rvect))*(hc/rc^4)*cy0+((2*hc)/(rc^2))*cdy0; 
    cddy = ((-mu/rc^3)+(hc^2/rc^4))*cy0 + 2*(dot(vvect,rvect))*(hc/rc^4)*cx0-2*(hc/rc^2)*cdx0;
    cddz = -(mu/rc^3)*cz0;
    % these are lineareized and have non circular assumption
    dstate_c = [cdx0; cdy0; cdz0; cddx; cddy; cddz];

    dstate = [dstate_t; dstate_c];

end