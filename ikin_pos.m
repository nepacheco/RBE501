function [t1,t2,t3,t4,t5,t6] = ikin_pos(desired_pos,desired_orient)
x = desired_pos(1);
y = desired_pos(2);
z = desired_pos(3);
% calculating the inverse position kinematics
t1 = atan2(y,x);

% solving for theta 2
d = sqrt(x^2 + y^2)-160;
r = sqrt(120^2 + 720^2);
a = sqrt(d^2 + (z-475)^2);
F = (-r^2 + a^2 + 600^2)/(2*a*600);
beta = atan2(sqrt(1-F^2),F);
gamma = atan2(z-475,d);
t2 = beta + gamma - pi/2;

% solving for theta 3
phi = atan2(120,720);
G = (-600^2 + a^2 + r^2)/(2*a*r);
delta = atan2(sqrt(1-G^2),G);
n = (pi/2 - gamma) + delta + phi;
t3 = pi/2-n-t2;

% calculating the inverse orientation kinematics
R04 =[-sin(t2 + t3)*cos(t1),  sin(t1), cos(t2 + t3)*cos(t1);
    -sin(t2 + t3)*sin(t1), -cos(t1), cos(t2 + t3)*sin(t1);
    cos(t2 + t3),        0,         sin(t2 + t3)];

R57 = (R04')*desired_orient;
% these are based on forward kinematics of a spherical wrist
t4 = atan2(R57(2,3),R57(1,3)); 
t5 = atan2(sqrt(1-R57(3,3)^2),R57(3,3));
t6 = atan2(R57(3,2),-R57(3,1));
end