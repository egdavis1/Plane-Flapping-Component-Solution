a2 = 15;
r1 = [-0.8, 0, 0];
r2 = [0.35*cosd(a2), 0.35*sind(a2), 0];

syms a3 a4;

eqn1 = -0.8 - 0.9*cosd(a4) + 0.35*cosd(15) + cosd(a3);
eqn2 = -0.9*sind(a4) + 0.35*sind(15) + sind(a3);

solution = solve([eqn1, eqn2], [a3, a4]);

a3Sol = solution.a3;
a4Sol = solution.a4;

r3 = [cosd(a3Sol(1)), sind(a3Sol(1)), 0];
r4 = [-0.9*cosd(a4Sol(1)), -0.9*sind(a4Sol(1)), 0];

w2 = 0:-10:-2000;

vAx = -w2*r2(2);
vAy = w2*r2(1);

w3 = (vAx - cotd(a4Sol(1) - 90)*vAy) / (cotd(a4Sol(1)-90)*r3(1) + r3(2));

vB = (vAy + w3*r3(1))/sind(a4Sol(1) - 90);

vBx = vB*cosd(a4Sol(1) - 90);
vBy = vB*sind(a4Sol(1) - 90);
vBVec = [vBx, vBy, 0];
w4 = vBx/r4(2);



%Graph w2 vs w3 and w4:
figure
plot(w2, w3, 'b')
hold on;

plot(w2, w4, 'r')
hold off;

grid on;
title('Graph of w2 vs w3 and w4');
xlabel('w2 [rad/sec]');
ylabel('w3, w4 [rad/sec]');
legend({'w3','w4'},'Location','northeast')


aAx = -(w2.*w2)*r2(1);
aAy = -(w2.*w2)*r2(2);

aBnx = (vB.*vB/norm(r4))*cosd(180-a4Sol(1));
aBny = (vB.*vB/norm(r4))*sind(180-a4Sol(1));

alpha3 = (aAx - w3.*w3*r3(1) - aAy*cotd(a4Sol(1)-90) + w3.*w3*r3(2)*cotd(a4Sol(1)-90) - aBny*cotd(a4Sol(1)-90) - aBnx) / ...
         (r3(1)*cotd(a4Sol(1) - 90) + r3(2));


aB = (aAy + alpha3*r3(1) - w3.*w3*r3(2) + aBny)/sind(a4Sol(1)-90);
alpha4 = -aB/norm(r4);

figure
plot(w2, alpha3)
hold on;
plot(w2, alpha4)
hold off;

grid on;
title('Graph of w2 vs alpha3 and alpha4');
xlabel('w2 [rad/sec]');
ylabel('alpha3, alpha4 [rad/sec^2]');
legend({'alpha3','alpha4'},'Location','northeast')
