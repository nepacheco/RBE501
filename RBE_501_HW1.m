%% Homework 1
clc; clear; close all;

%% Homogeneous Transformation Matrices
syms t1 t2 t3 t4 t5 t6
DH_table = [t1 330 75 -sym(pi)/2;
            t2-sym(pi)/2 0 300 0;
            t3 0 75 -pi/2;
            t4+sym(pi) 320 0 -sym(pi)/2;
            t5 0 0 sym(pi)/2;
            t6 80 0 0];
T01 = tdh(DH_table(1,:));
T12 = tdh(DH_table(2,:));
T23 = tdh(DH_table(3,:));
T34 = tdh(DH_table(4,:));
T45 = tdh(DH_table(5,:));
T56 = tdh(DH_table(6,:));
%% Composite Transformation

T0_6 = T01*T12*T23*T34*T45*T56;
T0_6 = simplify(T0_6,'Steps', 100);
pretty(T0_6)
input = deg2rad([0,75,30,135,-45,60]);
T0_6_val = subs(T0_6, [t1,t2,t3,t4,t5,t6],input);
T0_6_val = vpa(T0_6_val,4)

%% Plotting the robot

plot_robot([0,75,30,135,-45,60]);
title('Position of robot at given input')
xlabel('x');
ylabel('y');
zlabel('z');
snapnow
%% bonus

fanimator(@(t)animator(t));
title('Animated Plot Moving to Position')
xlabel('x');
ylabel('y');
zlabel('z');
view([45,45])
axis equal
grid on
playAnimation;
snapnow
%% Functions
function p = plot_robot(q)
% assumes q is in degrees
q = deg2rad(q);
t1 = q(1); t2 = q(2); t3 = q(3); t4 = q(4); t5 = q(5); t6 = q(6);
DH_table = [t1 330 75 -pi/2;
            t2-pi/2 0 300 0;
            t3 0 75 -pi/2;
            t4+pi 320 0 -pi/2;
            t5 0 0 pi/2;
            t6 80 0 0];
T01 = tdh(DH_table(1,:));
T12 = tdh(DH_table(2,:));
T23 = tdh(DH_table(3,:));
T34 = tdh(DH_table(4,:));
T45 = tdh(DH_table(5,:));
T56 = tdh(DH_table(6,:));
T02 = T01*T12;
T03 = T01*T12*T23;
T04 = T01*T12*T23*T34;
T05 = T01*T12*T23*T34*T45;
T0_tip = T01*T12*T23*T34*T45*T56;

x = [0 T01(1,4) T02(1,4) T03(1,4) T04(1,4) T05(1,4) T0_tip(1,4)];
y = [0 T01(2,4) T02(2,4) T03(2,4) T04(2,4) T05(2,4) T0_tip(2,4)];
z = [0 T01(3,4) T02(3,4) T03(3,4) T04(3,4) T05(3,4) T0_tip(3,4)];
p = plot3(x,y,z, '-o','MarkerFaceColor', 'b');
view([45,45])
axis equal
grid on
end

function p = animator(t)
initial_pos = zeros(1,6);
increment = [0,75,30,135,-45,60]./10;
q = initial_pos+increment*t;
p = plot_robot(q);
end
