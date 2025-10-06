function earthHandle = addEarth(radius, axHandle)
    % addEarth - Generates a 3D Earth with corrected orientation and adds it to a specified axes.
    %
    % Syntax:
    %   earthHandle = addEarth(radius, axHandle)
    %
    % Inputs:
    %   radius - Radius of the Earth (default: 6371 km).
    %   axHandle - Axes handle where Earth will be plotted (default: gca).
    %
    % Outputs:
    %   earthHandle - Handle to the Earth surface object.

    if nargin < 1
        radius = 6378; % Default Earth's radius in kilometers
    end
    if nargin < 2
        axHandle = gca; % Default to the current axes
    end

    % Load Earth texture
    earthTexture = imread('https://eoimages.gsfc.nasa.gov/images/imagerecords/57000/57730/land_ocean_ice_2048.png');
    earthTexture = flipud(earthTexture); % Flip texture vertically to correct poles

    % Create a sphere
    [lon, lat, Z] = sphere(200); % Higher resolution for smoother appearance

    % Scale the sphere to the specified radius
    X = radius * lon;
    Y = radius * lat;

    % Create surface object
    earthHandle = surf(axHandle, X, Y, radius * Z, 'EdgeColor', 'none', ...
        'FaceColor', 'texturemap', 'CData', earthTexture);

    % Adjust lighting and material properties
    lightangle(axHandle, -45, 30);
    lighting(axHandle, 'phong');
    material(axHandle, 'shiny');
    axis(axHandle, 'equal');
    %axis(axHandle, 'off'); % Optional: Hide axes for a cleaner look
end