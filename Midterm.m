%% Mid term for RBE 501
clc;clear;close all
%% 2) DH Table
syms t1 t2 t3 t4 t5 t6 L real
dh_table = [0 t1 0 0;
            t2 0 -L sym(pi)/2;
            t3+sym(pi)/2 0 L -sym(pi)/2;
            0 2*L+t4 -L sym(pi);
            t5+sym(pi)/2 0 L -sym(pi)/2;
            t6-sym(pi)/2 3*L 0 0];
dh_table_var = @(t1, t2, t3, t4, t5, t6)...
            [0 t1 0 0;
            t2 0 -100 pi/2;
            t3+pi/2 0 100 -pi/2;
            0 2*100+t4 -100 pi;
            t5+pi/2 0 100 -pi/2;
            t6-pi/2 3*100 0 0];

%% 3) Generating Homogeneous transformation matrix

T_total = get_fwdkin(dh_table,true);
T_tip = T_total(:,:,6);

%% 4) Calculating position of end effector in home position
Home = subs(T_tip,[t1 t2 t3 t4 t5 t6,L],[zeros(1,6),100])
dh_table_home = subs(dh_table,[t1 t2 t3 t4 t5 t6,L],[zeros(1,6),100]);
plot_robot(dh_table_home);
view(-45,-45)


%% 5) Showing a vector in end effector frame relative to body
ee_vector = [10;10;10];
ee_point_bframe = Home*[ee_vector;1]
% or 
ee_vector_bframe = Home(1:3,1:3)*ee_vector

%% 6) Inverse Kinematics
desired_point = [-350; 50; -350];
% Finding possible solutions if we fix joint 2 and joint 3
pos = subs(T_tip(1:3,4),[L t2 t3],[100 0 0])
eqn1 = pos(1) == desired_point(1);
eqn2 = pos(2) == desired_point(2);
eqn3 = pos(3) == desired_point(3);
solution_ikin = solve([eqn1 eqn2 eqn3],[t1 t4 t5]);
t1_vals = real(vpa(solution_ikin.t1))
t4_vals = real(vpa(solution_ikin.t4))
t5_vals = real(vpa(solution_ikin.t5))

result_pos1 = vpa(subs(pos,[t1 t4 t5],[t1_vals(1) t4_vals(1) t5_vals(1)]),4)
result_pos2 = vpa(subs(pos,[t1 t4 t5],[t1_vals(2) t4_vals(2) t5_vals(2)]),4)

% finding possible solutions without fixing joints
pos = simplify(subs(T_tip(1:3,4),[L t2],[100 0]),'Steps',40)
eqn1 = pos(1) == desired_point(1);
eqn2 = pos(2) == desired_point(2);
eqn3 = pos(3) == desired_point(3);
solution_ikin = solve([eqn1 eqn2 eqn3],[t1 t3 t4 t5],'ReturnConditions',true);
t1_vals = real(simplify(solution_ikin.t1,'Steps',20))
t3_vals = real(simplify(solution_ikin.t3,'Steps',20))
t4_vals = real(simplify(solution_ikin.t4,'Steps',20))
t5_vals = real(simplify(solution_ikin.t5,'Steps',20))
conditions = solution_ikin.conditions

% By nature, theta 6 does not affect the position of the end effector and
% so there will always be an infinite number of solutions to the inverse
% kinematics problem.
%% 7) Jacobian
z0 = [0;0;1]; p0 = [0;0;0];
z1 = T_total(1:3,3,1); p1 = T_total(1:3,4,1);
z2 = T_total(1:3,3,2); p2 = T_total(1:3,4,2);
z3 = T_total(1:3,3,3); p3 = T_total(1:3,4,3);
z4 = T_total(1:3,3,4); p4 = T_total(1:3,4,4);
z5 = T_total(1:3,3,5); p5 = T_total(1:3,4,5);
pe = T_total(1:3,4,6);
Jv = simplify([z0 cross(z1,pe-p1) cross(z2,pe-p2)...
    z3 cross(z4,pe-p4) cross(z5,pe-p5)],'Steps', 10);
Jw = simplify([zeros(3,1) z1 z2 zeros(3,1) z4 z5], 'Steps', 10);
J = simplify([Jv; Jw],'Steps',10)

%% Determining Singularities
% Only care about positional singularities
Jv = simplify(subs(Jv,L,100),'Steps',20)
% velocity jacobian
det_Jv = simplify(det(Jv*Jv'),'Steps',20);
% t2 does not affect singularity as it is not in the determinant of Jv*Jv'
% Therefore we will set t2 = 0
Jv1 = simplify(subs(Jv,t2,0));
%% find where x velcoties will be equal to 0
% end effector is point directly in line with joint 2 and 3
Jvx = Jv1(1,2:5) == zeros(1,4)
% therefore t3 == +- pi/2
th3 = -pi/2;
% therefore (100*cos(t5)-300sin(t5) == 0
th5 = atan(1/3);
% therefore (t4 + 200) == 0
th4 = -200;
rank_typ1_ex1 = rank(vpa(subs(Jv,[t1 t2 t3 t4 t5 t6],[0 0 th3 th4 th5 0]),4))
% alternate th3
th3 = pi/2
rank_typ1_ex2 = rank(vpa(subs(Jv,[t1 t2 t3 t4 t5 t6],[0 0 th3 th4 th5 0]),4))

dh_table_sing_ex1 = subs(dh_table,[t1 t2 t3 t4 t5 t6,L],[0 0 -pi/2 -200 th5 0,100]);
plot_robot(dh_table_sing_ex1);
%% find where y velocities will be equal to 0
% This means the end effector is in line with joint 2
Jvy = Jv1(2,2:5) == zeros(1,4)
% therefore 300cos(t5) + 100sin(t5) == 0
th5 = atan(-3);
% This means that -cos(t3)(t4 + 200) == 100 and there are infinite
% solutions to this equation
%one example
th3 = acos(-1/2.5); th4 = 50;
rank_typ2_ex1 = rank(vpa(subs(Jv,[t1 t2 t3 t4 t5 t6],[0 0 th3 th4 th5 0]),4))
% another example
th3 = 2*pi/3; th4 = 0;
rank_typ2_ex2 =rank(vpa(subs(Jv,[t1 t2 t3 t4 t5 t6],[0 0 th3 th4 th5 0]),4))

dh_table_sing_ex2 = subs(dh_table,[t1 t2 t3 t4 t5 t6,L],[0 0 2*pi/3 0 th5 0,100]);
plot_robot(dh_table_sing_ex2);
%% Inverse Velocity
tip_vel = [10;0;10];
Jv_val = subs(Jv,[L t1 t2 t3 t4 t5 t6],[100 zeros(1,6)])
joint_vel = vpa(pinv(Jv_val)*tip_vel,4)

rank_Jv_val = rank(Jv_val)
% Becauses the rank of Jv = 3 (full rank) at the home position and there
% are 5 columns, there are an infinite number of possible solutions to the 
% inverse velocity problem. The pinv solution finds an answer that 
% minimizes the norm of joint velocities, but is only one possible answer.
% By moving only joint 1, 3, and 4 for example, there are an infinite
% combination of joint velocities that would produce a tip velocity of
% [10;0;10].

disp("Need hand written stuff here");
%% functions
function T = get_fwdkin(dh_table,is_sym)
    rows = size(dh_table,1);
    if is_sym
        T = sym('T',[4,4,rows]);
    else
        T = zeros(4,4,rows);
    end
    for i = 1:rows
        if i == 1
            T(:,:,i) = tdh(dh_table(i,:));
        else
            T(:,:,i) = simplify(T(:,:,i-1)*tdh(dh_table(i,:)),'Steps',10);
        end
    end
end

function p = plot_robot(dh_table)
    T = get_fwdkin(dh_table,false);
    num_transforms = size(T,3);
    pos = [0;0;0];
    for i = 1:num_transforms
        pos = [pos T(1:3,4,i)];
    end
%     p = plot3(pos(3,:),-pos(2,:),pos(1,:),'Marker','o');
    p = plot3(pos(3,:),-pos(2,:),pos(1,:),'Marker','o','LineWidth',2);
    xlabel('z axis (mm)')
    ylabel('y axis (mm)')
    zlabel('x axis (mm)')
    axis equal
    grid on
end

function f = det_Jv_fun(t)
    t3 = t(1); t4 = t(2); t5 = t(3);
    d = @(t3, t4 , t5)...
        24000000*t4 + 24004200000000*sin(2*t5) + 9600000000*cos(t3) + 100000000*t4*cos(t3)^2 + 240000*t4^2*cos(t3) - 35999800*t4*cos(t3)^3 + 400*t4^3*cos(t3) - 43999600*t4*cos(t3)^4 + 128000000*t4*cos(t5)^2 + 48000000*t4*sin(2*t5) - 4800000000*cos(t5)*sin(t3) + 2399980000*sin(t3)*sin(t5) + 6200030000*cos(t3)^2 - 3999960000*cos(t3)^3 - 1999970000*cos(t3)^4 + 28011200000000*cos(t5)^2 - 28000000000000*cos(t5)^4 + 570000*t4^2*cos(t3)^2 - 120000*t4^2*cos(t3)^3 + 1600*t4^3*cos(t3)^2 - 269999*t4^2*cos(t3)^4 + 320000*t4^2*cos(t5)^2 - 200*t4^3*cos(t3)^3 + 2*t4^4*cos(t3)^2 - 800*t4^3*cos(t3)^4 - t4^4*cos(t3)^4 - 12800000000*cos(t3)*cos(t5)^2 + 120000*t4^2*sin(2*t5) - 96000000000000*cos(t5)^3*sin(t5) - 3999980000*sin(t3)^3*sin(t5) + 80000000*t4*cos(t3) - 19999840000*cos(t3)^2*cos(t5)^2 + 12800000000*cos(t3)^3*cos(t5)^2 + 5600000000*cos(t3)^2*cos(t5)^4 + 19999920000*cos(t3)^4*cos(t5)^2 - 5600000000*cos(t3)^4*cos(t5)^4 + 60000*t4^2 + 20000*t4^2*sin(t3)*sin(t5) - 24000000*t4*sin(t3)^3*sin(t5) - 2000000000*cos(t3)*sin(t3)*sin(t5) - 320000000*t4*cos(t3)^2*cos(t5)^2 + 64000000*t4*cos(t3)^3*cos(t5)^2 + 224000000*t4*cos(t3)^4*cos(t5)^2 + 9599940000*cos(t3)^2*cos(t5)*sin(t3) + 14400000000*cos(t3)*cos(t5)^3*sin(t3) + 10799880000*cos(t3)^3*cos(t5)*sin(t3) - 20399880000*cos(t3)^2*cos(t5)*sin(t5) + 9600000000*cos(t3)^3*cos(t5)*sin(t5) + 20399940000*cos(t3)^4*cos(t5)*sin(t5) - 60000*t4^2*sin(t3)^3*sin(t5) + 5999960000*cos(t3)^3*sin(t3)*sin(t5) - 48000000*t4*cos(t5)*sin(t3) - 800000*t4^2*cos(t3)^2*cos(t5)^2 + 560000*t4^2*cos(t3)^4*cos(t5)^2 + 8000000*t4*sin(t3)*sin(t5) - 7200000000*cos(t3)^2*cos(t5)^3*sin(t3) - 21600000000*cos(t3)^3*cos(t5)^3*sin(t3) + 19200000000*cos(t3)^2*cos(t5)^3*sin(t5) - 19200000000*cos(t3)^4*cos(t5)^3*sin(t5) - 64000000*t4*cos(t3)*cos(t5)^2 - 120000*t4^2*cos(t5)*sin(t3) - 1200000000*cos(t3)*cos(t5)*sin(t3) - 9600000000*cos(t3)*cos(t5)*sin(t5) + 180000*t4^2*cos(t3)^2*cos(t5)*sin(t3) + 720000*t4^2*cos(t3)^3*cos(t5)*sin(t3) - 108000000*t4*cos(t3)^3*cos(t5)^3*sin(t3) - 600000*t4^2*cos(t3)^2*cos(t5)*sin(t5) + 1200*t4^3*cos(t3)^3*cos(t5)*sin(t3) + 420000*t4^2*cos(t3)^4*cos(t5)*sin(t5) + 240000*t4^2*cos(t3)^3*sin(t3)*sin(t5) + 400*t4^3*cos(t3)^3*sin(t3)*sin(t5) + 20800000000*cos(t3)*cos(t5)^2*sin(t3)*sin(t5) - 102000000*t4*cos(t3)*cos(t5)*sin(t3) - 48000000*t4*cos(t3)*cos(t5)*sin(t5) - 42000000*t4*cos(t3)*sin(t3)*sin(t5) - 10400000000*cos(t3)^2*cos(t5)^2*sin(t3)*sin(t5) - 31200000000*cos(t3)^3*cos(t5)^2*sin(t3)*sin(t5) + 72000000*t4*cos(t3)^2*cos(t5)*sin(t3) - 720000*t4^2*cos(t3)*cos(t5)*sin(t3) + 72000000*t4*cos(t3)*cos(t5)^3*sin(t3) + 149999400*t4*cos(t3)^3*cos(t5)*sin(t3) - 1200*t4^3*cos(t3)*cos(t5)*sin(t3) - 240000000*t4*cos(t3)^2*cos(t5)*sin(t5) + 48000000*t4*cos(t3)^3*cos(t5)*sin(t5) + 168000000*t4*cos(t3)^4*cos(t5)*sin(t5) - 240000*t4^2*cos(t3)*sin(t3)*sin(t5) + 61999800*t4*cos(t3)^3*sin(t3)*sin(t5) - 400*t4^3*cos(t3)*sin(t3)*sin(t5) - 156000000*t4*cos(t3)^3*cos(t5)^2*sin(t3)*sin(t5) + 104000000*t4*cos(t3)*cos(t5)^2*sin(t3)*sin(t5) + 9004200000000;
    f = d(t3,t4,t5);
end