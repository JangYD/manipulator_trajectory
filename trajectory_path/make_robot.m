function make_robot = make_robot(dh_parameter)
    robot = rigidBodyTree;   
    robot_joint = size(dh_parameter, 1);
    bodies = cell(robot_joint, 1);
    joints = cell(robot_joint, 1);

    for i = 1:robot_joint
        bodies{i} = rigidBody(['body' num2str(i)]);
        joints{i} = rigidBodyJoint(['joint' num2str(i)], "revolute");
        setFixedTransform(joints{i}, dh_parameter(i,:),'dh');
        bodies{i}.Joint = joints{i};
        
        if i == 1
            addBody(robot, bodies{i}, "base");
        else
            addBody(robot, bodies{i}, bodies{i-1}.Name);
        end
    end
    make_robot = robot;
end