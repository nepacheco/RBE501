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