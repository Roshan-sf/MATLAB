function R = Cy(theta)
% R2d  Find the rotation matrix about the number two axis for a rotation of
% theta degrees
%

R = [cos(theta) 0 -sin(theta);
      0          1 0;
     sin(theta) 0 cos(theta)];