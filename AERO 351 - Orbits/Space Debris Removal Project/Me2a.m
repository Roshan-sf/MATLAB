function [a] = Me2a(Me,mu)
    %ME2A Finds Semi-Major Axis from Mean Motion
    %   [a] = Me2a(Me,mu)

    Me1 = Me*((2*pi)/86400); %Converts from rev/day to rad/s
    a = (mu/(Me1^2))^(1/3); %finds Semi-Major axis (a)
end