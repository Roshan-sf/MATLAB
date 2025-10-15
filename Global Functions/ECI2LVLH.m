function [Q, QT] = ECI2LVLH(R_ECI,V_ECI)
%Rotation Matrix for satellite relative positioning
%Usage: Q, QT = ECI2LVLH[R,V] where inputs are properties of chief/target
%Position should be a 3x1 col vector: Q*posVec_ECI = posVec_LVLH

    ha = cross(R_ECI,V_ECI);
    
    ihat = R_ECI/norm(R_ECI);
    khat = ha/norm(ha);
    jhat = cross(khat,ihat);

    Q = [ihat'; jhat'; khat'];
    QT = Q';

end