function tleData = readTLE(filename)
    % Reads a Two-Line Element (TLE) from a text file and parses the data.
    % 
    % Input:
    %   filename - Name of the text file containing the TLE
    % 
    % Output:
    %   tleData - A struct containing the extracted TLE information
    
    % Open the file for reading
    
    fileID = fopen(filename, 'r');
    if fileID == -1
        error('Error opening file: %s', filename);
    end
    
    % Read the two lines of the TLE
    line1 = fgetl(fileID);
    line2 = fgetl(fileID);
    fclose(fileID);
    
    % Check if lines are correctly read
    if length(line1) < 69 || length(line2) < 69
        error('Invalid TLE format');
    end
    
    % Parse line 1
    tleData.satelliteNumber = str2double(line1(3:7)); % Satellite Catalog Number
    tleData.classification = strtrim(line1(8));       % Classification
    tleData.intDesignator = strtrim(line1(10:17));    % International Designator
    tleData.epochYear = str2double(line1(19:20));     % Epoch Year
    tleData.epochDay = str2double(line1(21:32));      % Epoch Day
    tleData.meanMotionDerivative = str2double(line1(34:43)); % First Derivative of Mean Motion
    tleData.meanMotionSecDerivative = str2double(line1(45:50)) * 10^str2double(line1(51:52)); % Second Derivative
    tleData.bstar = str2double(line1(54:59)) * 10^str2double(line1(60:61)); % BSTAR Drag Term
    tleData.ephemerisType = str2double(line1(63));    % Ephemeris Type
    tleData.elementSetNumber = str2double(line1(65:68)); % Element Set Number
    
    % Parse line 2
    tleData.inclination = str2double(line2(9:16));    % Inclination (degrees)
    tleData.rightAscension = str2double(line2(18:25));% Right Ascension of Ascending Node (degrees)
    tleData.eccentricity = str2double(['0.', line2(27:33)]); % Eccentricity
    tleData.argumentOfPerigee = str2double(line2(35:42)); % Argument of Perigee (degrees)
    tleData.meanAnomaly = str2double(line2(44:51));   % Mean Anomaly (degrees)
    tleData.meanMotion = str2double(line2(53:63));    % Mean Motion (revolutions per day)
    tleData.revolutionNumber = str2double(line2(65:68)); % Revolution Number at Epoch
    
    % Validate year
    if tleData.epochYear < 57
        tleData.epochYear = tleData.epochYear + 2000; % 21st century
    else
        tleData.epochYear = tleData.epochYear + 1900; % 20th century
    end
end