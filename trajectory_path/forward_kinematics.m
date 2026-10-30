function [robot_end_position, link_position, R_matrix] = forward_kinematics(dh_parameter)
    n = size(dh_parameter, 1);
    robot_end_position = eye(4);
    link_position = cell(n, 1);
    R_matrix = cell(n, 1);

    for i = 1:n
        a = dh_parameter(i,1); alpha = dh_parameter(i,2);
        d = dh_parameter(i,3); theta = dh_parameter(i,4);

        cos_theta = cos(theta); 
        sin_theta = sin(theta);
        cos_alpha = cos(alpha); 
        sin_alpha = sin(alpha);

        A = [cos_theta, -sin_theta*cos_alpha,   sin_theta*sin_alpha,  a*cos_theta;
             sin_theta,  cos_theta*cos_alpha, -cos_theta*sin_alpha, a*sin_theta;
             0,                 sin_alpha,            cos_alpha,           d;
             0,                      0,                    0,                1];
        
        robot_end_position = robot_end_position * A;
        link_position{i} = robot_end_position;
        R_matrix{i} = A;
    end
end
