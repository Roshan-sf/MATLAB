    %% Roshan Jaiswal-Ferri
    %Section - 01 
    %Aero 355 Lab 5 Vibrations: 3/12/25
    
    %% Workspace Prep
    
    format long     %Allows for more accurate decimals
    close all;      %Clears all
    clear all;      %Clears Workspace
    clc;            %Clears Command Window
    
    %% Data Manip
    
    w_rod = 0.09;
    w_mass = 0.13;
    mass = (w_mass)/(32.2*12);
    massr = (w_rod)/(32.2*12);
    g = 32.2;
    d = 0.25; %in
    L = 11; %in
    x = 11;
    h = 11;
    w = w_rod/L;
    Irod = pi*(d)^4/64;
    E = 29000000; %psi
    P = w_mass;
    
    delta_rod = ((w*(x)^2)/(24*E*Irod)) * (x^2 + 6*L^2 - 4*L*x);
    delta_mass = ((P*x^2)/(6*E*Irod)) * (3*h - x);
    
    total_delta = delta_mass + delta_rod;
    F = w_rod + w_mass;
    
    %khard
    mass2 = mass + massr;
    %Q = 33.3;
    total_delta_fn = total_delta;%*Q;
    k_hard = F/total_delta_fn; %
    fn_hard = (1/(2*pi))*sqrt(k_hard/mass2); %
    
    %keasy
    k_easy = (3*E*Irod)/(L^3);
    fn_easy = (1/(2*pi))*sqrt(k_easy/mass);
    
    %Finding moment
    Q = 1/(2*0.015);
    M = Q*(((w*L^2)/2)+P*x);
    
    %Finding Stress
    c = d/2;
    sigma = (M*c)/Irod;
    
    %Displacement
    disph = Q*(delta_mass+delta_rod); %displacement hard
    dispe = Q*(delta_mass); %easy way
    
    %% Now everything again but empirical
    
    Data = readtable('2017-5-24_AERO353_vibelab.csv');
    accel = Data.Article13749_M_Filtered__g_(2:end); %in gs
    freq = Data.X_Data_Hz_(2:end);
    [maxa,pos] = max(accel);
    %[~,hz_at_max] = min(abs(accel-maxa));
    freqM = freq(pos);
    
    Q = max(accel/0.55);
    
    M2 = Q*(((w*L^2)/2)+P*x);
    
    %Finding Stress
    sigma2 = (M2*c)/Irod;
    
    %Displacement
    displ = Q*(delta_mass+delta_rod); %displacement empirical
    
    %% Error Analysis
    
    lsra = 0.000001/2;
    lsrx = 0.5;
    lsrl = 0.5; %use error of 0.5 for length
    lsrd = 0.01/2; 
    lsrw = 0.005+0.5;
    Eunc = 0.5;
    Iunc = lsrd^4;
    lsrp = lsrw;
    I = Irod;
    
    %%
    
    dtdx = (((3*w)/(24*E*I))*(4*x^3)) + (2*P*(3*x^2))/(6*E*I);
    
    dtdi = (-1/(I^2)) * (((w*x^2)/(24*E))*(3*x^2)+(((P*x^2)*2*x)/(6*E)));
    
    dtdw = (x^2/(24*E*I))*3*x^2;
    
    dtdp = x^3/(3*E*I);
    
    error_disp_hard = (dtdx*lsrx) + (dtdi*Iunc) + (dtdw*lsrw) + (dtdp*lsrp);
    
    %%
    
    F = P;
    
    dkdfh = 1/disph;
    dkddisph = -F/(2*disph^2);
    
    error_k_hard = (dkddisph*lsrw) + (dkdfh*error_disp_hard);
    
    %%
    
    dkdi = (3*E)/L^3; %Check this, it is huge
    dkdl = (-9*E*I)/L^4;
    
    error_k_easy = (dkdi*Iunc) + (dkdl*lsrl);
    
    %%
    
    for i = 1:2
        if i == 2
            Q = 33.3;
        end
    
        dmdq = ((w*L^2)/2) + (P*x);
        dmdx = P*Q + Q*w*x;
        dmdp = (Q*x);
        dmdw = (Q*x^2)/2;
        if i == 1
            error_moment_post = (dmdq*lsra) + (dmdx*lsrx) + (dmdp*lsrp) + (dmdw*lsrw); % Q???
        elseif i == 2
            error_moment_pre = (dmdq*0) + (dmdx*lsrx) + (dmdp*lsrp) + (dmdw*lsrw);
        end
        
    end
    %%
    
    dsdi = (-1*M*c)/(I^2);
    dsdm = c/I;                 %all bad :(
    dsdc = M/I;
    
    error_sigma_easy = (dsdi*Iunc) + (dsdm*error_moment_pre) + (dsdc*lsrd);
    error_sigma_hard = (dsdi*Iunc) + (dsdm*error_moment_post) + (dsdc*lsrd);
    
    %%
    
    dddx = (P*3*x^2) / (6*E*I);
    dddi = (-1/(I^2)) * ((P*x^2)/(6*E));
    dddp = x^3/(6*E*I);
    
    error_disp_easy = (dddx*lsrx) + (dddi*Iunc) + (dddp*lsrp);
    
    %%
    
    dfdk = (1/(4*pi))*(1/sqrt(k_hard*mass2));
    dfdm = ((-sqrt(k_hard))/(4*pi))*mass2^(-3/2);
    
    dfdke = (1/(4*pi))*(1/sqrt(k_easy*mass));
    dfdme = -(k_easy^(1/2)) / (4 * pi * mass^(3/2));
    
    error_freq_hard = (dfdk*error_k_hard) + (dfdm*lsrw);
    error_freq_easy = (dfdke*error_k_easy) + (dfdme*lsrw);
    
    
    
    x = 1;