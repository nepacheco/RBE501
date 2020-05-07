clc;clear; close all;
%% Trajectory Generation for single joint over 3 seconds
q = [20 70 100 120]*pi/180;
vel = [0 0 0 0];
accel = [0 0 0 0];
t = [0 1 2 3];

a = zeros(6,3);
for i = 1:3
    a(:,i) = traj_generator([q(i);vel(i);accel(i)],[q(i+1);vel(i+1);accel(i+1)],t(i),t(i+1));
end

joint_pos = [];
joint_vel = [];
joint_accel = [];
t_span = [];
for i = 1:3
    t = linspace(i-1,i,1000)';
    t_span = [t_span;t];
    t_vec = [ones(1000,1) t t.^2 t.^3 t.^4 t.^5];
    td_vec = [zeros(1000,1) ones(1000,1) 2*t 3*t.^2 4*t.^3 5*t.^4];
    tdd_vec = [zeros(1000,1) zeros(1000,1) 2*ones(1000,1) 6*t 12*t.^2 20*t.^3];
    
    joint_pos = [joint_pos;t_vec*a(:,i)];
    joint_vel = [joint_vel;td_vec*a(:,i)];
    joint_accel = [joint_accel;tdd_vec*a(:,i)];
end
figure()
plot(t_span,joint_pos)
title('Joint position over time');
xlabel("time (s)");
ylabel("Joint angle (rad)");

figure()
plot(t_span,joint_vel);
title('Joint velocity over time');
xlabel("time (s)");
ylabel("Joint velocity (rad/s)")

figure()
plot(t_span,joint_accel);
title('Joint acceleration over time');
xlabel("time (s)");
ylabel("Joint acceleration (rad/s^2)")
%% Two Link Arm Robot trajectory
l1 = 300; l2 = 300;
m1 = 0.5; m2 = 0.5;

x0 = [300;450];
q0 = ikin(x0);
xf = [-300;450];
qf = ikin(xf);
t = [0 5];
v = 0;

a1 = traj_generator([q0(1);v],[qf(1);v],t(1),t(2));
a2 = traj_generator([q0(2);v],[qf(2);v],t(1),t(2));

points = 1000;
t_span = linspace(t(1),t(2),points)';
t = t_span;
t_vec = [ones(points,1) t t.^2 t.^3];
td_vec = [zeros(points,1) ones(points,1) 2*t 3*t.^2];
tdd_vec = [zeros(points,1) zeros(points,1) 2*ones(points,1) 6*t];

q = zeros(points,2); q_dot = zeros(points,2); q_dotdot = zeros(points,2);
% getting position velocity and acceleration on joint data
q(:,1) = t_vec*a1; q(:,2) = t_vec*a2;
q_dot(:,1) = td_vec*a1; q_dot(:,2) = td_vec*a2;
q_dotdot(:,1) = tdd_vec*a1; q_dotdot(:,2) = tdd_vec*a2;
pos = fwkin_pos(q');

figure()
plot(t_span,q(:,1),'r-',t_span,q(:,2),'b-')
title('Joint position over time');
xlabel("time (s)");
ylabel("Joint angle (rad)");
legend('Joint1','Joint2');

figure()
plot(t_span,q_dot(:,1),'r-',t_span,q_dot(:,2),'b-');
title('Joint velocity over time');
xlabel("time (s)");
ylabel("Joint velocity (rad/s)")
legend('Joint1','Joint2');

figure()
plot(t_span,q_dotdot(:,1),t_span,q_dotdot(:,2));
title('Joint acceleration over time');
xlabel("time (s)");
ylabel("Joint acceleration (rad/s^2)")
legend('Joint1','Joint2');

figure()
plot(pos(1,:),pos(2,:));
title("Tip Position")
xlabel("X Axis (mm)");
ylabel("Z Axis (mm)");

%% Dynamical Model of the system
l1 = 300; l2 = 300;
m1 = 0.5; m2 = 0.5;
syms q1 q2 q1_dot q2_dot q1_dotdot q2_dotdot real

dh_table_sym = [q1 0 l1 0; q2 0 l2 0];

dh_table_m1 = dh_table_sym(1,:);
T0_m1_total = get_fwdkin(dh_table_m1,true);
T0_m1 = T0_m1_total;
Jv1 = [cross([0;0;1],T0_m1(1:3,4))';0 0 0; 0 0 0]';
Jv1 = Jv1(1:2,1:2);

% generating transform to mass 2
dh_table_m2 = dh_table_sym;
T0_m2_total = get_fwdkin(dh_table_m2,true);
T0_m2 = T0_m2_total(:,:,2);
% jacobian stuff
T0_1_m2 = T0_m2_total(:,:,1);
pe_m2 = T0_m2(1:3,4);
Jv2 = simplify([cross([0;0;1],pe_m2)';
       cross(T0_1_m2(1:3,3),pe_m2-T0_1_m2(1:3,4))';
       0 0 0]','Steps',10);
Jvtip = Jv2(1:3,1:2);
Jv2 = Jv2(1:2,1:2);
Jwtip  = [[0 0 1];T0_1_m2(1:3,3)']';
Jtip = [Jvtip;Jwtip];

Jv = sym('Jv',[2 2 2]);
Jv(:,:,1) = Jv1; Jv(:,:,2) = Jv2;

m = [0.5;0.5];
D = zeros(2,2);
for i = 1:2
    lin = simplify(m(i)*Jv(:,:,i)'*Jv(:,:,i),'Steps',30);
    D = simplify(D + lin,'Steps',20);
end

C = sym('C',[2 2]);
q = [q1;q2];
q_dot = [q1_dot;q2_dot];
for k = 1:2
    for j = 1:2
        C(k,j) = sym(0);
        for i = 1:2
            C(k,j) = C(k,j) + simplify(1/2*(diff(D(k,j),q(i)) +...
                diff(D(k,i),q(j)) -...
                diff(D(i,j),q(k)))*q_dot(i),'Steps',50);
        end
    end
end
g = 9.8;
g_vec = [0 g 0]';
P = 0;
pos = sym('pos',[3 2]);
pos(:,1) = T0_m1(1:3,4); pos(:,2) = T0_m2(1:3,4);
for i = 1:2
    P = P + m(i)*g_vec'*pos(:,i);
end
G = simplify([diff(P,q1); diff(P,q2)],'Steps',20);

B = sym('B',[2,2]);
Ftip = sym('Ftip',[6,1]);
q_dotdot = [q1_dotdot;q2_dotdot];
tau = D*q_dotdot + C*q_dot + G + B*q_dot + Jtip'*Ftip;
pretty(tau)

%% Block Diagram

%% Control Law
syms Kp1 Kp2 Kv1 Kv2 q1d q2d q1d_dot q2d_dot q1d_dotdot q2d_dotdot 
qd = [q1d;q2d];
qd_dot = [q1d_dot;q2d_dot];
qd_dotdot = [q1d_dotdot; q2d_dotdot];
Kp = [Kp1 0; 0 Kp2]; Kv = [Kv1 0 ; 0 Kv2];

tau_control = Kp*(qd-q) + Kv*(qd_dot - q_dot) + D*qd_dotdot + C*qd_dot + G + B*q_dot + Jtip'*Ftip;
pretty(tau_control)

disp("TODO FIX TAU_CONTROL TO USE QDESIRED IN MATRICES m,d,c")
%% Mobile Robot Constraints
syms a r theta real
syms wr wl x_dot y_dot theta_dot wf wsw real
zeta_dot = [x_dot;y_dot;theta_dot];
R = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];

alpha_r = -sym(pi)/2; beta_r = sym(pi);
R_right = [sin(alpha_r + beta_r) -cos(alpha_r + beta_r) (-a/2)*cos(beta_r)]*R*zeta_dot == r*wr
S_right = [cos(alpha_r + beta_r) sin(alpha_r + beta_r) a/2*sin(beta_r)]*R*zeta_dot == 0

alpha_l = sym(pi)/2; beta_l = 0;
R_left = [sin(alpha_l + beta_l) -cos(alpha_l + beta_l) (-a/2)*cos(beta_l)]*R*zeta_dot == r*wl
S_left = [cos(alpha_l + beta_l) sin(alpha_l + beta_l) a/2*sin(beta_l)]*R*zeta_dot == 0

alpha_f = sym(pi); beta_f = sym(pi/2); gamma_f = 0;
syms c
R_front = [sin(alpha_f + beta_f + gamma_f) -cos(alpha_f + beta_f + gamma_f) (c)*cos(beta_f + gamma_f)]*R*zeta_dot == r*wf
S_front = [cos(alpha_f + beta_f + gamma_f) sin(alpha_f + beta_f + gamma_f) c*sin(beta_f + gamma_f)]*R*zeta_dot == r*wf*sin(gamma_f) - r*wsw

%% Constraint Matrix
C_mat = [sin(alpha_r + beta_r) -cos(alpha_r + beta_r) (-a/2)*cos(beta_r);
    sin(alpha_l + beta_l) -cos(alpha_l + beta_l) (-a/2)*cos(beta_l);
    sin(alpha_f + beta_f + gamma_f) -cos(alpha_f + beta_f + gamma_f) (c)*cos(beta_f + gamma_f);
    cos(alpha_r + beta_r) sin(alpha_r + beta_r) a/2*sin(beta_r);
    cos(alpha_l + beta_l) sin(alpha_l + beta_l) a/2*sin(beta_l);
cos(alpha_f + beta_f + gamma_f) sin(alpha_f + beta_f + gamma_f) c*sin(beta_f + gamma_f)]*R*zeta_dot...
    == [r*wr;
    r*wl;
    r*wf;
    0;
    0;
    r*wf*sin(gamma_f) - r*wsw]

% Removing unessesary constraints
C_mat_opt = [C_mat(1:2);C_mat(5)]

%% Velocity Kinematics
syms r wl wr theta real
R = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
zetaR_dot = [r*wl/2 + r*wr/2; 0; -r*wl/a + r*wr/a];
zetaI_dot = R'*zetaR_dot

%% Combined Forward Kinematics
syms r xr yr theta real
syms px py pz real
R_RW = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
T_RT = [sym('R',[3,3],'real') [px;py;pz]; zeros(1,3) 1]
T_WR = [R_RW' [xr;yr;0];zeros(1,3) 1]

T_WT = T_WR*T_RT

%% Combined velocity Kinematics
J_TR = sym('J_TR',[6,3],'real');
syms q1 q2 q3 q1_dot q2_dot q3_dot real
q_dot = [q1_dot;q2_dot;q3_dot];
v_theta = [cross([0;0;1],T_WT(1:3,4)-T_WR(1:3,4));0; 0; 1]*zetaI_dot(3) 
total_tip_vel = v_theta + [zetaI_dot(1:2);0;0;0;0] + J_TR*q_dot

%% Combined Jacobian
J_total = [((r*cos(theta))/2 + (r*(py*cos(theta) + px*sin(theta)))/a), ((r*cos(theta))/2 - (r*(py*cos(theta) + px*sin(theta)))/a), J_TR(1,1:3);
((r*sin(theta))/2 - (r*(px*cos(theta) - py*sin(theta)))/a)  ((r*sin(theta))/2 + (r*(px*cos(theta) - py*sin(theta)))/a)  J_TR(2,1:3);
0 0 J_TR(3,1:3);
0 0 J_TR(4,1:3);
0 0 J_TR(5,1:3);
(-r/a) (r/a) J_TR(6,1:3)];

%% torque at each wheel calculation
syms fx fy fz nx ny nz real
F = [fx fy fz nx ny nz]';
tau = J_total'*F;
%% ikin test
q = ikin([300;450])
%% functions

function q = ikin(pos)
q = zeros(size(pos,1),size(pos,2));
for i = 1:size(pos,2)
    l1 = 300; l2 = 300;
    q = zeros(2,1);
    x = pos(1,i);
    z = pos(2,i);
    h = sqrt(x^2 + z^2);
    D = (l1^2 + l2^2 - h^2)/(2*l1*l2);
    gamma = atan2(sqrt(1-D^2),D);
    q(2,i) = pi - gamma;
    
    beta = atan2(z,x);
    E = (h^2 + l1^2 - l2^2)/(2*h*l1);
    alpha = atan2(sqrt(1-E^2),E);
    q(1,i) = beta - alpha;
end
end

function x = fwkin_pos(q)
l1 = 300; l2 = 300;
x = zeros(size(q,1),size(q,2));
for i = 1:size(q,2)
    x(1,i) = l1*cos(q(1,i)) + l2*cos(q(1,i)+q(2,i));
    x(2,i) = l1*sin(q(1,i)) + l2*sin(q(1,i)+q(2,i));
end
end

function T = tdh(dh)
theta = dh(1);
d = dh(2);
a = dh(3);
alpha = dh(4);
T = [cos(theta) -sin(theta)*cos(alpha) sin(theta)*sin(alpha) a*cos(theta);
    sin(theta) cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta);
    0           sin(alpha)             cos(alpha)              d;
    0              0                   0                       1];
end

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
        T(:,:,i) = T(:,:,i-1)*tdh(dh_table(i,:));
        if is_sym
            T(:,:,i) = simplify(T(:,:,i),'Steps',20);
        end
    end
end
end
