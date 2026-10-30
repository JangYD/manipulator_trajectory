function [q, q_dot, q_dot_dot] = cubic_polynomial(path, traj_time, via_time)
    [row,col] = size(path);
    n = col-1;

    A = [1          0               0                   0;
         0          1               0                   0;
         1     via_time(2)   via_time(2)^2    via_time(2)^3;
         0          1           2*via_time(2)    3*via_time(2)^2];

    path_v = spline(path,via_time);
    path_v = [zeros(row,1) path_v zeros(row,1)];

    coeffs = zeros(row,4,n);
    fun = zeros(row,4,n);
    for i = 1 : n
        time = via_time(i+1) - via_time(i);

        a1 = path(:,i);
        a2 = path_v(:,i);
        a3 = ( 3*(path(:,i+1) - path(:,i)) - (path_v(:,i+1) +2*path_v(:,i))*time )/time^2;
        a4 = ( 2*(path(:,i) -path(:,i+1)) + (path_v(:,i+1) + path_v(:,i))*time )/time^3;

        for j = 1 : row
            coeffs(j,:,i) = [a1(j); a2(j); a3(j); a4(j)];
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
            via_q{j,i} = @(t) coeffs_temp(1) + coeffs_temp(2)*t + coeffs_temp(3)*t.^2 + coeffs_temp(4)*t.^3;
            via_q_dot{j,i} = @(t) coeffs_temp(2) + 2*coeffs_temp(3)*t + 3*coeffs_temp(4)*t.^2;
            via_q_dot_dot{j,i} = @(t) 2*coeffs_temp(3) + 6*coeffs_temp(4)*t;
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
        idx = loop_ind(i,loop_ind(i,:) > 0);
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
    v = zeros(path_row,waypoint_length);
    
    %%%%%%%%%%%%%%%%%%%%%doolittle_LUP_decomposition%%%%%%%%%%%%%%%%%
    % b = zeros(path_row,waypoint_length);
    % y = v;
    % T_length = size(b,2);
    % T = eye(T_length)*4;
    % side_T =ones(T_length-1,1);
    % T = T + diag(side_T,1) + diag(side_T,-1);
    % [L,U,P] = doolittle_LUP_decomposition(T);
    % 
    % 
    % for i = 1 : waypoint_length
    %     path_time = time(i+1) -time(i);
    %     b(:,i) = (3/path_time)*(path(:,i+2) - path(:,i));
    % end 
    % 
    % for i = 1 : size(b,1)
    %     temp = P*b(i,:)';
    %     b(i,:) = temp;
    % end
    % 
    % for i = 1 : size(b,1)
    %     temp = forward_substitution(L,b(i,:));
    %     y(i,:) = temp';
    % 
    % end
    % 
    % for i = 1 : size(b,1)
    %     temp = back_substitution(U,y(i,:));
    %     v(i,:) = temp';
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%% thomas %%%%%%%%%%%%%%%%%%%%%%%%
    RHS = zeros(path_row,path_col-2);
    [RHS_row,RhS_col] = size(RHS);
    T_length = RhS_col;
    T = eye(T_length)*4;
    side_T =ones(T_length-1,1);
    T = T + diag(side_T,1) + diag(side_T,-1);
    RHS = zeros(path_row,waypoint_length);
    for i = 1 : path_col-2
        traj_time = time(i+1)-time(i);
        RHS(:,i) = (3/traj_time)*(path(:,i+2) - path(:,i));
    end
    for i = 1 : RHS_row
        [~,~,x]=Thomas_Algorithm(T,RHS(i,:));
        v(i,:) = x;
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


% function y = forward_substitution(L_matrix,b)
%     [row, col] = size(b');
%     A = L_matrix;
%     y = zeros(row,col);
%     for i = 1 : row
%         if i == 1 
%             y(i) = b(i);
%             A(i+1:end,i) = A(i+1:end,i)*y(i);
%             continue
%         end
%         y(i) = b(i) - sum(A(i,1:i-1));
%         if i == row
%             break
%         end
%         A(i+1:end,i) = A(i+1:end,i)*y(i);
% 
%     end
% 
% end
% 
% function x = back_substitution(U_matrix,y)
%     [row, col] = size(y');
%     A = U_matrix;
%     x = zeros(row,col);
% 
%     for i = row : -1 : 1
%         if i == row
%             x(i) = y(i)/A(i,i);
%             A(1:i-1,i) =  A(1:i-1,i)*x(i);
%             continue
%         end
%         x(i) = (y(i) - sum(A(i,i+1:end)))/A(i,i);
%         if i == 1
%             break
%         end
%         A(1:i-1,i) =  A(1:i-1,i)*x(i);
%     end
% end
