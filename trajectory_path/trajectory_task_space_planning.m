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
sampling_time = 0.05;
  
% waypoints = [0 0 0;
%             0.15 0.2 0.1; 
%             -0.1 -0.1 0.1;  
%             -0.2 -0.3 0.2 ; 
%             0 0 0]'+robot_end_matrix(1:3,4);
% waypoints = [
%     % -0.25   -0.25   -0.15
%     -0.10   -0.15   -0.05;
%      0.00   0.00   0.00;
%      0.10   0.15   0.05;
%      0.20   0.20   0.10;
%      ]';
% waypoints = [...
%     0.00   0.00   0.00;
%     0.10   0.12   0.03;
%     0.20   0.20   0.10;
%     0.30   0.12   0.20;
%     0.40   0.00   0.30]' + robot_end_matrix(1:3,4);

waypoints = [ ...
   -0.2933   -0.1433   -0.3933   -0.4233   -0.2933;
   -0.2907   -0.0907   -0.3907   -0.4907   -0.2907;
    0.6696    0.7696    0.7696    0.5696    0.6696
]+0.08;

section_time = 4;
via_point_time = 0:section_time:(size(waypoints,2)-1)*section_time;

via_point_vel = zeros(size(waypoints));

traj_time = 0:sampling_time:via_point_time(end);
via_pos =  waypoints;
orient(1).target_quat=[0 -0.7071 0 -0.7071];
orient(1).tolerance = deg2rad(0.1);
  
%%%%%%%%%%%%%%  trajectory   %%%%%%%%%%%%%%%%%%%%
trajectory_type = 's_curve';
switch trajectory_type 
    case 'cubic'
        [q, dq, ddq] = cubic_polynomial(via_pos, traj_time, via_point_time);
    case 'quintic'
        [q, dq, ddq] = quintic_polynomial(via_pos, traj_time, via_point_time);
    case 'trapezoidal'
        [q, dq, ddq] = trapezoidal_velocity(via_pos, traj_time, via_point_time);
    case 's_curve'
        [q, dq, ddq] = s_curve(via_pos, traj_time, via_point_time);
end

%%%%%%%%%%%%%%  graph   %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
video_path = fullfile(pwd, 'task_space_quintic.mp4');

video = VideoWriter(video_path, 'MPEG-4');
video.FrameRate = 30;
video.Quality = 100;
open(video);

%%%%%%%%%%%%%%%%%%%%%%%%
start_q = th;
figure; hold on;
axis equal;
grid on
view([50 38])
p_traj=plot3(via_pos(1,1),via_pos(2,1),via_pos(3,1),'b.-');
plot3(via_pos(1,:),via_pos(2,:),via_pos(3,:),'ro--');
target_pos = q;
for i = 1 : numel(traj_time)
    t_start = tic;
    new_theta = LM_opt_inverse_kinematics(dh_param, target_pos(:,i), start_q, orient);
    elapsed_time = toc(t_start);
    start_q = new_theta;
    fprintf('Step %d / %d: IK 계산 시간 = %.6f초\n', i, numel(traj_time), elapsed_time);
    set(p_traj, 'xdata', [p_traj.XData target_pos(1,i)], ...
                'ydata', [p_traj.YData target_pos(2,i)], ...
                'zdata', [p_traj.ZData target_pos(3,i)]);
    
    show(robot,start_q,'Frames','on','PreservePlot',false);
    title(['Time :' num2str(traj_time(i)) 'sec']);
    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');
    drawnow
    
    check_dh_param =dh_param;
    check_dh_param(:,4) = start_q';
    end_pos = forward_kinematics(check_dh_param);
    %%%%%%%%%%%%%%%%%%%%%%%%
    frame = getframe(gcf);
    writeVideo(video, frame);
    %%%%%%%%%%%%%%%%%%%%%%%%
end

close(video);

figure(3);
subplot(3,1,1);
plot(traj_time,q);
xlabel('t(sec)');
ylabel('q(t) (rad)');
title('linear segments with parabolic blends (lspb)')
grid on
subplot(3,1,2);
plot(traj_time, dq);
xlabel('t (sec)');
ylabel('dq/dt(t) (rad/s)');
grid on
subplot(3,1,3);
plot(traj_time, ddq);
xlabel('t (sec)');
ylabel('d^{2}q/dt^{2}(t) (rad/s^{2})')
grid on

