function [q, q_dot, q_dot_dot] = trapezoidal_velocity(path, traj_time, via_time)
    [row,col] = size(path);
    n = col-1;
    
    q = zeros(row, size(traj_time,2));
    q_dot = zeros(row, size(traj_time,2));
    q_dot_dot = zeros(row, size(traj_time,2));

    via_section_time = (via_time(end) - via_time(1))/n;
    dt = traj_time(2)- traj_time(1);
    diff_q = diff(path,1,2);
    diff_t = diff(via_time);
    tb = round((diff_t/3)/dt)*dt;
    v = diff_q./(diff_t -tb);
    tb = tb(1);

    loop_ind = zeros(3,1);
    for i = 1 : 3
        if i == 1
            len = find((via_time(i) <= traj_time) & (traj_time <= tb));
            loop_ind(i,1:length(len)) = len;
        elseif i == 3
            len = find((via_section_time - tb < traj_time) & (traj_time <= via_section_time));
            loop_ind(i,1:length(len)) = len;
        else
            len = find((tb < traj_time) & (traj_time <= via_section_time - tb));
            loop_ind(i,1:length(len)) = len;
        end
    end
    
    time_loop_ind =  zeros(n,1);
    for i = 1 : n
        if i == 1
            len = find((via_time(i) <= traj_time) & (traj_time <= via_time(i+1)));
            time_loop_ind(i,1:length(len)) = len;
        else
            len = find((via_time(i) < traj_time) & (traj_time <= via_time(i+1)));
            time_loop_ind(i,1:length(len)) = len;
        end
    end

    for i = 1 : n
        idx = time_loop_ind(i, time_loop_ind(i,:) > 0);
        tau = traj_time(idx) - via_time(i);

        tau1 = find((0 <= tau) & (tau <= tb));
        tau2 = find((tb < tau) & (tau <= via_section_time - tb));
        tau3 = find((via_section_time - tb < tau) & (tau <= via_section_time));

        q1 = path(1:end,i) + (v(1:end,i).*tau(tau1).^2)/(2*tb);
        q2 = path(1:end,i) - v(1:end,i)*tb/2 + v(1:end,i)*tau(tau2);
        q3 = path(1:end,i+1) - v(1:end,i)*tau(tau3).^2/(2*tb) - v(1:end,i)*(via_section_time).^2/(2*tb) + v(1:end,i)*tau(tau3)*(via_section_time)/tb;

        q1_dot = (v(1:end,i).*tau(tau1))/tb;
        q2_dot = repmat(v(1:end,i),1,size(tau2,2));
        q3_dot = - v(1:end,i)*tau(tau3)/tb + (v(1:end,i)*(via_section_time))/tb;

        q1_dot_dot = repmat((v(1:end,i))/tb, 1,size(tau1,2));
        q2_dot_dot = repmat(zeros(row,1),1,size(tau2,2));
        q3_dot_dot = repmat(- v(1:end,i)/tb, 1,size(tau3,2));

        q(:,idx) = [q1 q2 q3];
        q_dot(:,idx) = [q1_dot q2_dot q3_dot];
        q_dot_dot(:,idx) = [q1_dot_dot q2_dot_dot q3_dot_dot];
    end 
end