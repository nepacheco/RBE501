%%3001 test script
DH_table = [t1 100 0 sym(pi)/2;
    t2 0 100 0;
    t3 0 120 0];
T01 = tdh(DH_table(1,:));
T12 = tdh(DH_table(2,:));
T23 = tdh(DH_table(3,:)); 
T02 = T01*T12;
T03 = T02*T23;
pe = T03(1:3,4);
Jv = [cross([0;0;1],pe) cross(T01(1:3,3),pe-T01(1:3,4)) cross(T02(1:3,3),pe-T02(1:3,4))];
% geometric jacobian
Jw = [[0;0;1] T01(1:3,3) T02(1:3,3)];

Jv = simplify(Jv,'Steps',10);
detJv = simplify(det(Jv),'Steps',10)

subs(detJv,[t1,t2,t3],[0,pi/12,0])