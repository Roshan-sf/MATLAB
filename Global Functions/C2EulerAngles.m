function euler_angles = C2EulerAngles(C)
% Converts a rotation matrix C to 3-2-1 Euler angles [phi; theta; psi]
% where:
% - phi = roll (rotation about x)
% - theta = pitch (rotation about y)
% - psi = yaw (rotation about z)

    theta = -asin(C(3,1));  % pitch
    if abs(cos(theta)) > 1e-6  % avoid gimbal lock
        phi = atan2(C(3,2), C(3,3));    % roll
        psi = atan2(C(2,1), C(1,1));    % yaw
    else
        % Gimbal lock case
        phi = 0;
        psi = atan2(-C(1,2), C(2,2));
    end

    euler_angles = [phi; theta; psi];
end
