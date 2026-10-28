clc;
clear;
%%%%%%%%%%% robot data %%%%%%%%%%%%%%%%
a = [0 -0.6127 -0.57155 0 0 0];
al = [pi/2 0 0 pi/2 -pi/2 0];
d = [0.1807 0 0 0.17415 0.11985 0.11655];
th = [0 0 0 0 0 0];
inv_q = sym('q', [1 6], 'real');

dh_param = dh_parameter(a, al, d, th);


robot = make_robot(dh_param);
robot.DataFormat = 'row';

[robot_end_matrix, robot_link_matrix] = forward_kinematics(dh_param);
robot_length = length(robot_link_matrix);

%%%%%%%%%%%%% Robot Start Pos  %%%%%%%%%%%%%%%%%%%%
for i = 1 : robot_length
    robot.Bodies{i}.Joint.HomePosition = th(i);   
end

%%%%%%%%%%%%%%  Create Via  %%%%%%%%%%%%%%%%%%%%
sampling_time = 0.1;
% waypoints = [0 0 0;
%             -0.25 0.2 0.1; 
%             -0.2 -0.1 0.3;  
%             -0.3 -0.3 0.2 ; 
%             0 0 0]';
% waypoints = [
%     -0.25   -0.25   -0.15
%     -0.10   -0.15   -0.05;
%      0.00   0.00   0.00;
%      0.10   0.15   0.05;
%      % 0.25   0.25   0.15;
%      ]';

waypoints = [ ...
   -0.2933   -0.1433   -0.3933   -0.4933   -0.2933;
   -0.2907   -0.0907   -0.3907   -0.5907   -0.2907;
    0.6696    0.7696    0.7696    0.8696    0.6696
];
section_time = 4;
via_point_time = 0:section_time:(size(waypoints,2)-1)*section_time;

via_point_vel = zeros(size(waypoints));

traj_time = 0:sampling_time:via_point_time(end);
via_pos =  waypoints;
orient(1).target_quat=[0 -0.7071 0 -0.7071];
orient(1).tolerance = deg2rad(0.1);

%%%%%%%%%%%%%%  inverse kinematics  %%%%%%%%%%%%%%%%%%%%
joint_value = zeros(robot_length,size(via_pos,2));
start_q = th;
target_pos = via_pos;
for i = 1 : size(via_pos,2)
    new_theta = LM_opt_inverse_kinematics(dh_param, target_pos(:,i), start_q, orient);
    start_q = new_theta;
    joint_value(:,i) = new_theta';
end

%%%%%%%%%%%%%%  trajectory   %%%%%%%%%%%%%%%%%%%%
trajectory_type = 'quintic';
switch trajectory_type 
    case 'cubic'
        [q, dq, ddq] = cubic_polynomial(joint_value, traj_time, via_point_time);
    case 'quintic'
        [q, dq, ddq] = quintic_polynomial(joint_value, traj_time, via_point_time);
    case 'trapezoidal'
        [q, dq, ddq] = trapezoidal_velocity(joint_value, traj_time, via_point_time);
    case 's_curve'
        [q, dq, ddq] = s_curve(joint_value, traj_time, via_point_time);
end

%%%%%%%%%%%%%%  graph   %%%%%%%%%%%%%%%%%%%%

figure; hold on;
axis equal;
grid on
view([50 38])

p_traj=plot3(via_pos(1,1),via_pos(2,1),via_pos(3,1),'b.-');

plot3(via_pos(1,:),via_pos(2,:),via_pos(3,:),'ro--');

for i = 1 : numel(traj_time)
    t_start = tic;
    forward_dh_para = dh_parameter(a, al, d, q(:,i)');
    elapsed_time = toc(t_start);
    fprintf('Step %d / %d: IK 계산 시간 = %.6f초\n', i, numel(traj_time), elapsed_time);

    end_pos = forward_kinematics(forward_dh_para);
    set(p_traj, 'xdata', [p_traj.XData end_pos(1,4)], ...
                'ydata', [p_traj.YData end_pos(2,4)], ...
                'zdata', [p_traj.ZData end_pos(3,4)]);
    % plotTransforms(end_pos(1:3,4)',repmat([1 0 0 0],[size(q(:,i),2) 1]),'FrameSize',0.03);
    show(robot,q(:,i)','Frames','on','PreservePlot',false);
    title('Trajectory Waypoints');
    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');
    drawnow
end


figure(3);
subplot(3,1,1);
plot(traj_time,q);
xlabel('t(sec)');
ylabel('q(t) (rad)');
title('linear segments with parabolic blends (lspb)')
subplot(3,1,2);
plot(traj_time(1:end), dq);
xlabel('t (sec)');
ylabel('dq/dt(t) (rad/s)');
subplot(3,1,3);
plot(traj_time(1:end), ddq);
xlabel('t (sec)');
ylabel('d^{2}q/dt^{2}(t) (rad/s^{2})')



