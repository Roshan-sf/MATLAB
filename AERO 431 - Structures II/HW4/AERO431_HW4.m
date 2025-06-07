%% Alessandro Tedeschi, Stefan Rosu, & Roshan Jaiswal-Ferri
%Section - 01
%Aero 431 HW4: 6/1/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Question 1

% Constants
E = 70e9;             % Young modulus [Pa]
G = 26e9;             % Shear modulus [Pa]
sig_y = 280e6;        % Yield strength [Pa]
nu = (E / (2 * G)) - 1; % Poisson ratio

t = 2e-3;            % Thickness [m]
T = 5e-3;            % Flange thickness [m]
h = 0.1;              % Web height [m]
w = 0.05;             % Flange width [m]

% Initiailize L and b 
L_ = linspace(0.1, 0.5, 6);
b_ = linspace(0.1, 0.3, 6);

% Simply supported plate buckling coefficient 
k_plate = @(AR, m) (pi^2) * (m^2 + (AR * m)^2);

% Simply supported plate buckling stress
buckling_stress = @(E, k, t_, b_) (k * pi^2 * E / (12 * (1 - nu^2))) * (t_ / b_)^2;

% Prepare result array
results = [];

% For loop to iterate through all L and b possibilities
for i = 1:length(L_)

for j = 1:length(b_)

    L = L_(i);
    b = b_(j);
    mode = 1;

    % Calculate aspect ratios for each structural element
    AR_skin = L / b;
    AR_web = h / b;
    AR_flange = L / (w / 2);
    AR_Iplate = h / b;

    % Calculate respective buckling coefficients from aspect ratio
    k_skin = k_plate(AR_skin, mode);
    k_web = k_plate(AR_web, mode);
    k_flange = k_plate(AR_flange, mode);
    k_Iplate = k_plate(AR_Iplate, mode);

    % Calculate buckling stress for each element
    sigma_skin = min(buckling_stress(E, k_skin, t, b), sig_y);
    sigma_web = min(buckling_stress(E, k_web, t, b), sig_y);
    sigma_flange = min(buckling_stress(E, k_flange, T, w/2), sig_y);
    sigma_Iplate = min(buckling_stress(E, k_Iplate, T, b), sig_y);
        I_beam = 2 * (w * T^3) / 12 + 2 * w * T * (h/2)^2;
        A_beam = (2 * w * T) + (h * t);
    sigma_Euler = min((pi^2 * E * I_beam) / (A_beam * L^2), sig_y);

    % Minimum stresses as critical stress value 
    sigmas_plate = [sigma_skin, sigma_web, sigma_flange, sigma_Iplate];
    sigmas_Ibeam = [sigma_skin, sigma_web, sigma_flange, sigma_Euler];

    [sigma_crit_plate, i_plate] = min(sigmas_plate);
    [sigma_crit_Ibeam, i_Ibeam] = min(sigmas_Ibeam);

    components = {'Skin panel', 'Spar web', 'Spar flange', 'I-beam plate'};

    if i_Ibeam == 4
        components{4} = 'I-beam';
    end

    % Store results
    results = [results;
        L, b, ...
        sigma_skin, sigma_web, sigma_flange, sigma_Iplate, sigma_Euler, ...
        sigma_crit_plate, sigma_crit_Ibeam, ...
        i_plate, i_Ibeam];
end

end

% Display as table
headers = {'L [m]', 'b [m]', 'sig_skin [Pa]', 'sig_web [Pa]', 'sig_flange [Pa]', ...
           'sig_Iplate [Pa]', 'sig_Euler [Pa]', 'sig_crit_plate [Pa]', ...
           'sig_crit_Ibeam [Pa]', 'Mode_Plate_Idx', 'Mode_Ibeam_Idx'};

% Organize results
solution = array2table(results, 'VariableNames', headers);
mode_names = {'Skin Panel', 'Spar Web', 'Spar Flange', 'I-beam Plate/Beam'};
solution.Mode_Plate = mode_names(solution.Mode_Plate_Idx)';
solution.Mode_Ibeam = mode_names(solution.Mode_Ibeam_Idx)';

% Display results
disp(solution)

% The mode numbers seem to vary slightly, & the buckling stress seems to 
% decreases with rib spacing. It was also oberved that the panel was the 
% first to buckle. Increasing rib / spar count will improve strength.

%% Problem 2

% Part 1
syms a % crack length [m]

w = 0.5; % specimen width [m]
sigma = 50e6; % applied stress [Pa]

% Define alpha = a/W
alpha = a / w;

% Define the beta function (geometry factor)
beta = (1.122 - 1.122*alpha - 0.820*alpha^2 + 3.768*alpha^3 - 3.040*alpha^4) / sqrt(1 - 2*alpha);

% Stress intensity factor K_I
K1 = beta * sigma * sqrt(pi * a); % [Pa·sqrt(m)]

K1c = 24e6; % MPa*sqrt(m) 

eqn = K1 == K1c;

a = min(real(double(solve(eqn, a))));
disp("Critical Crack Length: " + a +" m")

% Part 2
m = 3.59;
C = 3.15e-11;
sig_max = 50e6;       % Maximum stress [Pa]

a0 = a/2;             % Initial crack length [m]
ac = a0*2;

dsig = 50; % MPa (keep these units for emirical relation)
dN = 100;

Ncycles = 0; % Initialize
Ncycles2 = 0;
a = a0;
while a < ac
    alpha1 = a / w;
    beta1 = (1.122 - 1.122*alpha1- 0.820*alpha1^2 + 3.768*alpha1^3 - 3.040*alpha1^4) / sqrt(1 - 2*alpha1);
    dK = beta1*dsig*sqrt(pi*a);

    a = a + C*dN*dK^m;
    Ncycles = Ncycles + dN;
end

disp("Number of Cycles(dN = 100): " + Ncycles)

dN = 10;
a = a0;
while a < ac
    alpha1 = a / w;
    beta1 = (1.122 - 1.122*alpha1- 0.820*alpha1^2 + 3.768*alpha1^3 - 3.040*alpha1^4) / sqrt(1 - 2*alpha1);
    dK = beta1*dsig*sqrt(pi*a);
    a = a + C*dN*dK^m;
    Ncycles2 = Ncycles2 + dN;
end

disp("Number of Cycles(dN = 10): " + Ncycles2)

% Comment:
% The difference between the dN = 100 and dN = 10 cycles is not very large.
% (only about 80 cycles, which is negligible)

%% Question 3

% i) With a paper with a central crack, the size of the crack seems to 
% dictate the critical stress value that the paper can hold before the 
% crack propagates. Additionally, the sheet of paper fails by the crack 
% growing from the center of the paper towards the edges. We found that 
% the tear was very uniform and continued in the direction of the crack 
% nearly perfectly.
% 
% ii) The folded paper with flanges/stringer design and a central crack 
% appears to be more resistant to crack propagation and failure from 
% tensile stress. The crack of the folded paper took more tensile force 
% to grow than the flat sheet of paper. This means that the folded paper 
% has a greater critical tensile stress than the flat paper. The crack 
% still propagates from the inside to the outside similar to the flat 
% paper, and maintains the same direction as the original central crack. 
% The tear appears slightly more straight and uniform than the first trial 
% with the flat sheet of paper, this is most likely due to the increase in 
% rigidity from the bends in the paper or “stringers.” 
% 
% iii) We conducted the experiment with each configuration one time:
% 
% •	Configuration 1: narrow glue surface, ~1 cm wide 
% •	Configuration 2: wide glue surface, ~3 cm wide
% 
% We found that the wider glue surface caused the tensile load to be 
% better transferred between each sheet of paper and increase the overall 
% resistance to the crack propagation. 
% 
% •	Configuration 3: weak glue
% •	Configuration 4: strong glue
% 
% Having a higher quality, or stronger, glue had more resistance to crack 
% propagation. This is a similar result as having a wider glue surface. 
% We suspect that the strength of the bond between sheets of paper is 
% directly related to the resistance against tearing. The bond strength 
% is increased whether the strength of the glue is increased or the 
% surface area of the bond is increased. 
% 
% We found that for each configuration, after the center sheet of paper 
% completely ripped, the glue would then fail and the outside sheets of 
% paper would not tear.
% 
% iv) We conducted the experiment with each configuration one time:
% 
% •	Configuration 1: pins close to the edge, ~1 cm
% •	Configuration 2: pins far from the edge, ~3 cm
% 
% We noticed that having pins closer to the edge of the paper causes the 
% edge of the paper to rip near the fasteners at smaller tensile loads. 
% 
% •	Configuration 3: less frequent pin placement
% •	Configuration 4: more frequent pin placement
% 
% Having more pins along the fastening surface increased the resistance to 
% crack propagation and increased the critical stress of the paper around 
% the fasteners. More pins means the load can be distributed to more 
% fasteners which lowers the individual stress on each one, meaning the 
% system as a whole can take a higher load.
% 
