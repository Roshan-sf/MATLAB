function [avgx] = avg(x,l)
%AVG Summary of this function goes here
%   Detailed explanation goes here
g = 0;

for i = 1:l
    g = g + x(1,i);
end

avgx = g/l;