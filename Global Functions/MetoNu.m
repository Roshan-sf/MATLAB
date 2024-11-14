function [theta] = MetoNu(Me, ecc)
%TA2TsincePeri Finds true anomaly when given time since perigee passage
%   [theta] = MetoNu(Me, ecc)
%   ecc = eccentricity (0<ecc<1)
%   t = time since perigee passage (sec)
%   T = Period (sec)
%   Me = 2*pi*t/T; % rads

% disp("Me: " + Me)
% disp("Heart: Me > pi!")

% Newtons method
if Me < pi
    Eold = Me + ecc/2;
else
    Eold = Me - ecc/2;
end

Enew = Eold - (Me-Eold+ecc*sin(Eold))/(-1+ecc*cos(Eold));
iter = 0;
while iter < 10000 & abs(Eold-Enew)>10^-8
    Eold = Enew;
    Enew = Eold - (Me-Eold+ecc*sin(Eold))/(-1+ecc*cos(Eold));
    iter = iter + 1;
end
theta = 2*atan(tan(Enew/2)*sqrt((1+ecc)/(1-ecc)));
end