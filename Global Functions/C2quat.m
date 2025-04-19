function q = C2quat(C)
% Converts a 3x3 rotation matrix to a quaternion [q1; q2; q3; q4]
% Output is a 4x1 quaternion, with q4 being the scalar part

    tr = trace(C);
    q = zeros(4,1);

    if tr > 0
        S = sqrt(tr + 1.0) * 2; % S = 4*q4
        q(4) = 0.25 * S;
        q(1) = (C(3,2) - C(2,3)) / S;
        q(2) = (C(1,3) - C(3,1)) / S;
        q(3) = (C(2,1) - C(1,2)) / S;
    elseif (C(1,1) > C(2,2)) && (C(1,1) > C(3,3))
        S = sqrt(1.0 + C(1,1) - C(2,2) - C(3,3)) * 2; % S = 4*q1
        q(4) = (C(3,2) - C(2,3)) / S;
        q(1) = 0.25 * S;
        q(2) = (C(1,2) + C(2,1)) / S;
        q(3) = (C(1,3) + C(3,1)) / S;
    elseif C(2,2) > C(3,3)
        S = sqrt(1.0 + C(2,2) - C(1,1) - C(3,3)) * 2; % S = 4*q2
        q(4) = (C(1,3) - C(3,1)) / S;
        q(1) = (C(1,2) + C(2,1)) / S;
        q(2) = 0.25 * S;
        q(3) = (C(2,3) + C(3,2)) / S;
    else
        S = sqrt(1.0 + C(3,3) - C(1,1) - C(2,2)) * 2; % S = 4*q3
        q(4) = (C(2,1) - C(1,2)) / S;
        q(1) = (C(1,3) + C(3,1)) / S;
        q(2) = (C(2,3) + C(3,2)) / S;
        q(3) = 0.25 * S;
    end
end
