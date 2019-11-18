w2 = -46;
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

syms vB w3;

w2Vec = [0, 0, w2];

vA = cross(w2Vec, r2);
w3Vec = [0, 0, w3];
vBVec = [vB*cosd(a4Sol(1)-90), vB*sind(a4Sol(1)-90), 0];

eqnVB = vA + cross(w3Vec, r3)- vBVec;

solutionVB = solve(eqnVB, [vB, w3]);

w3Sol = solutionVB.w3;
vBSol = solutionVB.vB;
double(vBSol);
double(w3Sol);

vBVec = [vBSol*cosd(a4Sol(1)-90), vBSol*sind(a4Sol(1)- 90), 0];
vBVec_new = [double(vBSol*cosd(a4Sol(1)-90)), double(vBSol*sind(a4Sol(1)- 90)), 0];
w3Vec = [0, 0, w3Sol];

syms w4;
w4Vec = [0, 0, w4];
eqnW4 = cross(w4Vec, (-1)*r4) - vBVec;
solutionW4 = solve(eqnW4(1), w4);

w4Vec = [0, 0, solutionW4];

%Solve for aA components
syms aAx aAy;

aAVec = [aAx, aAy, 0];

eqnAA = - w2*w2*r2 - aAVec;
solutionAA = solve(eqnAA, [aAx, aAy]);
aAxSol = solutionAA.aAx;
aAySol = solutionAA.aAy;

aAVec = [aAxSol, aAySol, 0];
double(aAVec);

syms aB alpha3;

alpha3Vec = [0, 0, alpha3];

aBVec = [aB*cosd(a4Sol(1) - 90) + (norm(vBVec)*norm(vBVec)/norm(r4))*cosd(180-a4Sol(1)), aB*sind(a4Sol(1) - 90) - (norm(vBVec)*norm(vBVec)/norm(r4))*sind(180-a4Sol(1)), 0];
eqnAB = aAVec + cross(alpha3Vec, r3) - norm(w3Vec)*norm(w3Vec)*r3 - aBVec;

solutionAB = solve(eqnAB, [aB, alpha3]);
aBSol = solutionAB.aB;
alpha3Sol = solutionAB.alpha3;

aBVec = [aBSol*cosd(90 - a4Sol(1)) + norm(vBVec)*norm(vBVec)/norm(r3)*cosd(180-a4Sol(1)), aBSol*sind(90 - a4Sol(1)) - norm(vBVec)*norm(vBVec)/norm(r3)*sind(180-a4Sol(1)), 0];
alpha3Vec = [0, 0, alpha3Sol];
double(aBVec);
double(alpha3Vec);

alpha4 = -aBSol/norm(r4);

fprintf('Theta 3 = %.3f\n', a3Sol(1));
fprintf('Theta 4 = %.3f\n', a4Sol(1));
fprintf('\n');
fprintf('vA = %.3f, %.3f, %.3f\n', vA(1), vA(2), vA(3));
fprintf('vB = %.3f, %.3f, %.3f\n', vBVec(1), vBVec(2), vBVec(3));
fprintf('\n');
fprintf('w3 = %.3f\n', w3Sol);
fprintf('w4 = %.3f\n', solutionW4);
fprintf('\n');
fprintf('aA = %.3f, %.3f, %.3f\n', aAVec(1), aAVec(2), aAVec(3));
fprintf('aB = %.3f, %.3f, %.3f\n', aBVec(1), aBVec(2), aBVec(3));
fprintf('\n');
fprintf('Alpha 3 = %.3f\n', alpha3Sol);
fprintf('Alpha 4 = %.3f\n', alpha4);
