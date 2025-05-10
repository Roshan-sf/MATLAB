% function pos = lla2eci(lat, lon, alt, JD)
%     % LLA2ECI  Convert geodetic latitude, longitude, altitude, and Julian date
%     % to Earth-Centered Inertial (ECI) coordinates.
%     %
%     % Inputs:
%     %   lat - latitude in degrees
%     %   lon - longitude in degrees
%     %   alt - altitude in kilometers
%     %   JD  - Julian date
%     %
%     % Output:
%     %   pos - 3x1 vector in ECI frame (km)
% 
%     % WGS-84 ellipsoid parameters
%     Re = 6378.137;              % Equatorial radius (km)
%     f = 1 / 298.257223563;      % Flattening
%     e2 = f * (2 - f);           % Square of eccentricity
% 
%     % Convert degrees to radians
%     lat_rad = deg2rad(lat);
%     lon_rad = deg2rad(lon);
% 
%     % Prime vertical radius of curvature
%     N = Re / sqrt(1 - e2 * sin(lat_rad)^2);
% 
%     % ECEF coordinates
%     x_ecef = (N + alt) * cos(lat_rad) * cos(lon_rad);
%     y_ecef = (N + alt) * cos(lat_rad) * sin(lon_rad);
%     z_ecef = (N * (1 - e2) + alt) * sin(lat_rad);
% 
%     % Greenwich Sidereal Time (GST) in degrees
%     gst = ct2lst(0, JD);  % Always use 0 to get GST
% 
%     % Rotate ECEF to ECI using GST
%     theta = deg2rad(gst);
%     R = [cos(theta), -sin(theta), 0;
%          sin(theta),  cos(theta), 0;
%          0,           0,          1];
% 
%     % ECI position
%     pos = R * [x_ecef; y_ecef; z_ecef];
% end

function pos = lla2eci(lat, lon, alt, JD)

Re = 6378.137;
gst = ct2lst(lon, JD); %deg
theta = (gst) * pi/180.0; %rad

r = (Re + alt)*cos(lat*pi/180.0); % km

pos(1,1) = r*cos(theta);% km
pos(2,1) = r*sin(theta);% km
pos(3,1) = (alt+Re)*sin(lat*pi/180.0);% km