%Roshan Jaiswal-Ferri
%Aero 215

function [SideC] = LawOfCosines(SideA, SideB,intAngAB)
    %Inputs are:SideA, SideB, and the interior angle
    % IMPORTANT: The interior angle is expected to be given in degrees not
    % radians
    %Output is the length of the unknown side

    %Roshan Jaiwal-Ferri
    %Using the Law of Cosines/the actual math:
    
    SideC = sqrt((SideA^2+SideB^2)-2*SideB*SideA*cosd(intAngAB));
end