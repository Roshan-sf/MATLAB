function [J0, JD] = tojd(day,month,year,ut)

%Floor is the integer portion of x, Y is the 4-digit year, M is the 2-digit month, and D
%is the 2-digit day of the month

timeV = sscanf(ut, '%d:%d:%d');
hours = timeV(1);
minutes = timeV(2);
seconds = timeV(3);

deciTime = hours + (minutes/60) + (seconds/3600);


J0 = 367*year - floor((7*(year+floor((month+9)/12)))/4) +floor((275*month)/9) + day + 1721013.5;
JD = J0 + (deciTime/24);

end