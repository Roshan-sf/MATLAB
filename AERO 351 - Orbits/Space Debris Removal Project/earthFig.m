% Load Earth texture
earthTexture = imread('https://eoimages.gsfc.nasa.gov/images/imagerecords/57000/57730/land_ocean_ice_2048.png');

% Create a sphere
[lon, lat, Z] = sphere(4); % Increase resolution by changing the number

% Scale the sphere to Earth's approximate dimensions
radius = 6378; % Earth's mean radius in kilometers
X = radius * lon; 
Y = radius * lat;

% Display the sphere with texture
figure;
earth = surf(X, Y, radius * Z, 'EdgeColor', 'none'); % Create 3D surface
hold on;

% Map texture to the sphere
set(earth, 'FaceColor', 'texturemap', 'CData', earthTexture);

% Set axis and lighting
axis equal;
%axis off;
lightangle(-45, 30); % Light source position
lighting phong;
material shiny;

% Add title
title('3D Model of Earth');
