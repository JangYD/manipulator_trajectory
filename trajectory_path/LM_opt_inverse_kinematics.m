function q = LM_opt_inverse_kinematics(dh_parameter, target_pos, q_guess, orient)
    max_iter = 1500;
    lambda = 1e-3;  
    nu = 10;
    curr_x = q_guess(:);
    grad_tol = 1e-10;
    sol_tol = 1e-10;
    pos_err_tol = 1e-10;
    if isempty(orient.tolerance)
        ori_err_tol = deg2rad(0.1);
    else
        ori_err_tol = orient.tolerance;
    end

    for i = 1 : max_iter
        j = jacobian_approximation(@(q) residual(dh_parameter, q, target_pos,  orient.target_quat), curr_x);
        r = residual(dh_parameter, curr_x, target_pos,  orient.target_quat);
        err_pos = r(1:3);
        err_ori = r(4:6);
        if norm(err_pos) < pos_err_tol && norm(err_ori) < ori_err_tol
            break;
        end
        % H = j'*j + lambda*eye(size(j));               % levenberg algorithm
        H = j'*j + lambda*diag(diag(j'*j));         % levenberg marquardt algorithm

        g = j'*r;
        if norm(g) < grad_tol
            break;
        end

        next_x = curr_x - H\g;
        if norm(next_x - curr_x) < sol_tol
            break;
        end

        next_r = residual(dh_parameter, next_x, target_pos,  orient.target_quat);

        if norm(next_r)^2 < norm(r)^2
            lambda = lambda/nu;
            curr_x = next_x;
        else
            lambda = lambda*nu;
        end        
        
    end
    curr_x = curr_x + 2*pi*round((q_guess(:) - curr_x)/(2*pi));
    q = curr_x';
end

function J = jacobian_approximation(f,x)
    eps = 1e-6;
    x = x(:);
    fx = f(x);
    col = length(x);
    row = length(fx);
    J = zeros(row,col);

    for i = 1 : col
        dx = zeros(col,1);
        dx(i) = eps;
        fx1 = f(x + dx);   % f(x + dx) - f(x)
        fx2 = f(x - dx);   % f(x - dx) - f(x)
        J(:,i) = (fx1 - fx2)/(2*eps);
    end
end

function error = residual(dh_parameter, q, target_pos, target_quat)
    q = q(:);
    J = jacobian(dh_parameter, q);
    pos = J.end_matrix(1:3,4);
    rot = J.end_matrix(1:3,1:3);

    pos_err = pos - target_pos;
    rot_err = quat2rotm(target_quat)'*rot;
    axis_angle_error = angle_axis_representation(rot_err);
    orient_error = axis_angle_error(1:3)*axis_angle_error(4);

    error = [pos_err(:); orient_error(:)];
end

function axis_angle = angle_axis_representation(rotation_matrix)
    axis_angle = zeros(1,4);
    theta = double(acos(1/2*(rotation_matrix(1,1) + rotation_matrix(2,2) + rotation_matrix(3,3) -1)));
    if abs(sin(theta)) > 1e-12
        axis_angle(1:3) = (1/(2*sin(theta)))*(vee(rotation_matrix - rotation_matrix'));
        axis_angle(end) = theta;
    else
       [v,D] = eig(eye(3) - rotation_matrix);
       axis_angle(1,1:3) = v(:,end)';
       axis_angle(end) = theta;   
    end
end

function w = vee(matrix)
    w = [matrix(3,2); matrix(1,3); matrix(2,1)];
end

function J = jacobian(dh_parameter, q_guess)   
    length = size(dh_parameter, 1);
    matrix = eye(4);
    pos = zeros(3,length+1);
    rot = zeros(3,length+1);
    pos(:,1) = [0;0;0];
    rot(:,1) = [0;0;1];
    for i = 1 : length
        q_guess_theta = q_guess(i,1);
        a = dh_parameter(i,1);
        al = dh_parameter(i,2);
        d = dh_parameter(i,3);

        A = Rz(q_guess_theta)*Tz(d)*Tx(a)*Rx(al);
        matrix = matrix*A;

        pos(:,i+1) = matrix(1:3,4);
        rot(:,i+1) = matrix(1:3,3);
    end
    J.pos = zeros(3,length);
    J.rot = zeros(3,length);

    for i = 1 : length
        J.pos(:,i) = pos(:,i+1);
        J.rot(:,i) = rot(:,i);
    end
    J.matrix(1:3 ,1:length) = J.pos;
    J.matrix(4:6 ,1:length) = J.rot;
    J.end_matrix = matrix;
end