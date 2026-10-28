function [q, q_dot, q_dot_dot] = quintic_polynomial(path, traj_time, via_time)
    syms t 'real';
    [row,col] = size(path);
    n = col-1;

    A = [1          0               0                   0           0              0;            
         0          1               0                   0           0              0;
         0          0               2                   0           0              0;
         1     via_time(2)   via_time(2)^2    via_time(2)^3   via_time(2)^4     via_time(2)^5;  
         0          1        2*via_time(2)   3*via_time(2)^2  4*via_time(2)^3   5*via_time(2)^4
         0          0               2          6*via_time(2)  12*via_time(2)^2  20*via_time(2)^3];
    
    path_v_tmp = spline(path,via_time);
    path_v = zeros(row,size(path_v_tmp,2)/2);
    path_a = zeros(row,size(path_v_tmp,2)/2);
    for i = 1 : size(path_v_tmp,2)/2
        path_v(1:row,i) = path_v_tmp(:,i*2-1);
        path_a(1:row,i) = path_v_tmp(:,i*2);
    end
    path_v = [zeros(row,1) path_v zeros(row,1)];
    path_a = [zeros(row,1) path_a zeros(row,1)];

    coeffs = zeros(row,6,n);
    fun = zeros(row,6,n);
    for i = 1 : n
        time = via_time(i+1) - via_time(i);
        
        a1 = path(:,i);
        a2 = path_v(:,i);
        a3 = 0.5*path_a(:,i);
        a4 =(20*(path(:,i+1) - path(:,i)) - (8*path_v(:,i+1) +12*path_v(:,i))*time + (path_a(:,i+1) - 3*path_a(:,i))*time^2 )/(2*time^3);
        a5 = (30*(path(:,i) - path(:,i+1)) +(14*path_v(:,i+1) +16*path_v(:,i))*time + (3*path_a(:,i) -2*path_a(:,i+1))*time^2 )/(2*time^4);
        a6 = (12*(path(:,i+1)- path(:,i)) - (6*path_v(:,i+1) + 6*path_v(:,i))*time + (path_a(:,i+1) - path_a(:,i))*time^2 )/(2*time^5);
        
        for j = 1 : row
            coeffs(j,:,i) = [a1(j); a2(j); a3(j); a4(j); a5(j); a6(j)];
        end

        for j = 1 : row
            fun_temp = A*coeffs(j,:,i)';
            fun(j,:,i) = fun_temp';
        end
        
    end
    
    via_q = cell(row,n);
    via_q_dot = cell(row,n);
    via_q_dot_dot = cell(row,n);

    for i = 1 : n
        for j = 1 : row
            coeffs_temp = coeffs(j,:,i);
            via_q{j,i} = @(t) coeffs_temp(1) + coeffs_temp(2)*t + coeffs_temp(3)*t.^2 + coeffs_temp(4)*t.^3 + coeffs_temp(5)*t.^4 + coeffs_temp(6)*t.^5;
            via_q_dot{j,i} = @(t) coeffs_temp(2) + 2*coeffs_temp(3)*t + 3*coeffs_temp(4)*t.^2 + 4*coeffs_temp(5)*t.^3 + 5*coeffs_temp(6)*t.^4;
            via_q_dot_dot{j,i} = @(t) 2*coeffs_temp(3) + 6*coeffs_temp(4)*t + 12*coeffs_temp(5)*t.^2 + 20*coeffs_temp(6)*t.^3;
        end
    end

    loop_ind = zeros(n,1);
    for i = 1 : n
        if i == 1
            len = find((via_time(i) <= traj_time) & (traj_time <= via_time(i+1)));
            loop_ind(i,1:length(len)) = len;
        else
            len = find(via_time(i) < traj_time & traj_time <= via_time(i+1));
            loop_ind(i,1:length(len)) = len;
        end
    end
        
    q = zeros(row, size(traj_time,2));
    q_dot = zeros(row,size(traj_time,2));
    q_dot_dot = zeros(row,size(traj_time,2));

    for i = 1 : n
        idx = loop_ind(i, loop_ind(i,:) > 0);
        tau = traj_time(idx) - via_time(i);
        
        for j = 1 : row
            q(j,idx) = via_q{j,i}(tau);
            q_dot(j,idx) = via_q_dot{j,i}(tau);
            q_dot_dot(j,idx) = via_q_dot_dot{j,i}(tau);
        end
    end    
end

function v = spline(matrix,time)
    
    path = matrix;
    [path_row,path_col] = size(path);
    waypoint_length = path_col-2;
    v = zeros(3,waypoint_length*2);
    y = v;
    Ajd = ones(waypoint_length,1);
    Ajtd = ones(waypoint_length-1,1);
    Asd = ones(waypoint_length,1);
    Astd = ones(waypoint_length-1,1);

    Aj = -4*(diag(Ajtd,-1)) + 4*(diag(Ajtd,1));
    Bj =  3*(diag(Ajd)) -(1/2)*(diag(Ajtd,-1)) - (1/2)*(diag(Ajtd,1));
    As = -16*diag(Asd) - 7*diag(Astd,-1) -7*diag(Astd,1);
    Bs =  zeros(waypoint_length,waypoint_length) - diag(Astd,-1) +diag(Astd,1);
    t = time(2)- time(1);
    A_tmp = [Aj*2 t*Bj*2; As t*Bs];
    A = A_tmp;
    for i = 1 : size(A_tmp,1)/2
        A(:,i+(i-1)) = A_tmp(:,i);
        A(:,i+i) = A_tmp(:,i+size(A_tmp,1)/2);
    end
    
    [L,U,P] = doolittle_LUP_decomposition(A);
    b =zeros(path_row,waypoint_length*2);
    
    [L_row,L_col] = size(L);
    for i = 1 : L_row/2
        b(:,i) = (20/(time(i+1)-time(i)))*(path(:,i+2) + path(:,i) - 2*path(:,i+1));
    end
    
    for i = 1 : L_col/2
        b(:,i+3) = (15/(time(i+1)-time(i)))*(path(:,i)- path(:,i+2));
    end
    
    for i = 1 : size(b,1)
        temp = P*b(i,:)';
        b(i,:) = temp';
    end
    
    for i = 1 : size(b,1)
        temp = forward_substitution(L,b(i,:));
        y(i,:) = temp';
        
    end
    
    for i = 1 : size(b,1)
        temp = back_substitution(U,y(i,:));
        v(i,:) = temp';
    end
end

function y = forward_substitution(L_matrix,b)
    [row, col] = size(b');
    A = L_matrix;
    y = zeros(row,col);
    for i = 1 : row
        if i == 1 
            y(i) = b(i);
            A(i+1:end,i) = A(i+1:end,i)*y(i);
            continue
        end
        y(i) = b(i) - sum(A(i,1:i-1));
        if i == row
            break
        end
        A(i+1:end,i) = A(i+1:end,i)*y(i);

    end

end

function x = back_substitution(U_matrix,y)
    [row, col] = size(y');
    A = U_matrix;
    x = zeros(row,col);

    for i = row : -1 : 1
        if i == row
            x(i) = y(i)/A(i,i);
            A(1:i-1,i) =  A(1:i-1,i)*x(i);
            continue
        end
        x(i) = (y(i) - sum(A(i,i+1:end)))/A(i,i);
        if i == 1
            break
        end
        A(1:i-1,i) =  A(1:i-1,i)*x(i);
    end
end

