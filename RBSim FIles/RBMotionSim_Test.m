% RB Sim Test Case

close all
clear

psi = 0:.1:4*pi;
theta = zeros(1, length(psi));
phi = zeros(1, length(psi));

RBMotionSim(psi, theta, phi);

close all
RBMotionSim(theta, psi, phi);

close all
RBMotionSim(phi, theta, psi);