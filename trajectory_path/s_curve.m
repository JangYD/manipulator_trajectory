function [q, qd, qdd] = s_curve(path,traj_time,via_time)
    [row,col] = size(path);
    qf = zeros(row,col-1);
    vmax = zeros(row,col-1);
    amax = zeros(row,col-1);
    jmax = zeros(row,col-1);
    delta_t = zeros(1,col-1);
    traj_T = zeros(1,col-1);
    T = zeros(1,col-1);
    

    for i = 1 : col-1
        qf(:,i) = path(:,i+1) - path(:,i);
        traj_T(i) = via_time(i+1) - via_time(i);
        dt = traj_time(2) - traj_time(1);
        delta_t(i) = round(traj_T(i)/7*(1/dt))*dt;
        vmax(:,i) = qf(:,i)/(traj_T(i) -2*delta_t(i) - delta_t(i));
        amax(:,i) = vmax(:,i)/(delta_t(i) + delta_t(i));
        jmax(:,i) = amax(:,i)/delta_t(i);
        T(i) = amax(1,i)./jmax(1,i) +vmax(1,i)./amax(1,i) +qf(1,i)./vmax(1,i);
    end
   
    t = zeros(7,1);
    t(1) = amax(1,1)/jmax(1,1);
    t(3) = t(1) + vmax(1,1)/amax(1,1);
    t(2) = t(3) - t(1);
    t(4) = T(i) - t(3);
    t(5) = T(i) - t(2);
    t(6) = T(i) - t(1);
    t(7) = t(6) + t(1);

    
    
    a = zeros(row,7,col-1);
    b = zeros(row,7,col-1);
    c = zeros(row,7,col-1);
    d = zeros(row,7,col-1);
    for i = 1 : col-1
        a(:,1,i) = jmax(:,i);
        a(:,2,i) = 0;
        a(:,3,i) = -jmax(:,i);
        a(:,4,i) = 0;
        a(:,5,i) = -jmax(:,i);
        a(:,6,i) = 0;
        a(:,7,i) = jmax(:,i);

        b(:,1,i) = 0;
        b(:,2,i) = amax(:,i);
        b(:,3,i) = amax(:,i) + jmax(:,i)*t(2);
        b(:,4,i) = 0;
        b(:,5,i) = jmax(:,i)*t(4);
        b(:,6,i) = -amax(:,i);
        b(:,7,i) = -amax(:,i) - jmax(:,i)*t(6);
    end
   
    c(:,1) = 0;
    for j = 1 : col-1
        d(:,1,j) = path(:,j);
        for i = 2 : 7
            c(:,i,j) = (1/2)*(a(:,i-1,j) - a(:,i,j))*t(i-1)^2 + (b(:,i-1,j) - b(:,i,j))*t(i-1) + c(:,i-1,j);
            d(:,i,j) = (1/6)*(a(:,i-1,j) - a(:,i,j))*t(i-1)^3  + (1/2)*(b(:,i-1,j) - b(:,i,j))*t(i-1)^2 + (c(:,i-1,j)-c(:,i,j))*t(i-1) + d(:,i-1,j);
        end
    end

    sampling_t = traj_time;
    seg_ind = zeros(7,1);
    
    for i = 1 : 7
        if i == 1
            len =find((0 <= sampling_t) & (sampling_t <= t(i)));
            seg_ind(i,1:length(len)) = len;
        else
            len = find((t(i-1) < sampling_t) & (sampling_t <= t(i)));
            seg_ind(i,1:length(len)) = len;
        end
    end

    q_fun = cell(7,col-1);
    qd_fun = cell(7,col-1);
    qdd_fun = cell(7,col-1);
    for j = 1 : col-1
        for i = 1 : 7
            q_fun{i,j} = @(tau) 1/6*a(:,i,j)*tau.^3 + 1/2*b(:,i,j)*tau.^2 +c(:,i,j)*tau +d(:,i,j);
            qd_fun{i,j} = @(tau) 1/2*a(:,i,j)*tau.^2 + b(:,i,j)*tau + c(:,i,j);
            qdd_fun{i,j} = @(tau) a(:,i,j)*tau + b(:,i,j);
        end
    end

    q = zeros(row,length(sampling_t));
    qd = zeros(row,length(sampling_t));
    qdd = zeros(row,length(sampling_t));

    for j = 1 : col-1
        for i = 1 : 7
            if j == 1
                idx = seg_ind(i, seg_ind(i,:) > 0);
            else
                seg_ind_min = min(seg_ind(1,1));
                seg_ind_max = max(seg_ind(7,:));
                cycle_value = seg_ind_max - seg_ind_min;
                q_cycle_first_idx = seg_ind(i,1) + cycle_value*(j-1);
                q_cycle_end_idx = max(seg_ind(i,:)) + cycle_value*(j-1);
                idx = q_cycle_first_idx:q_cycle_end_idx;

            end
            tau = traj_time(idx) - via_time(j);
            q(:,idx) = q_fun{i,j}(tau);
            qd(:,idx) = qd_fun{i,j}(tau);
            qdd(:,idx) = qdd_fun{i,j}(tau);
        end
    end
end