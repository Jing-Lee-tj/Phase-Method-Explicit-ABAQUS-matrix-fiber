syms C11 C22 C12 C23 C44;
syms E1 E2 V12 V23 G12;
C = [C11,C12,C12,0,0,0;
    C12,C22,C23,0,0,0;
    C12,C23,C22,0,0,0;
    0,0,0,C44,0,0;
    0,0,0,0,1/2*(C22-C23),0;
    0,0,0,0,0,C44];
CINV = C^-1;
% 定义方程组
eq1 = 1/E1 - CINV(1,1);
eq2 = -V12/E1 - CINV(1,2);
eq3 = 1/E2 - CINV(2,2);
eq4 = -V23/E2 - CINV(2,3);
eq5 = G12 -C44;
% 求解方程组
[C11_sol, C22_sol, C12_sol, C23_sol,C44_sol] = solve([eq1, eq2, eq3, eq4,eq5], [C11, C22, C12, C23,C44]);
% 验证
C_sol = [C11_sol,C12_sol,C12_sol,0,0,0;
    C12_sol,C22_sol,C23_sol,0,0,0;
    C12_sol,C23_sol,C22_sol,0,0,0;
    0,0,0,C44_sol,0,0;
    0,0,0,0,1/2*(C22_sol-C23_sol),0;
    0,0,0,0,0,C44_sol];
E1 = 26.5e9; %Longitudinal Young%s modulus
E2 = 2.6e9;  %Transverse Young’s modulus
G12 = 1.3e9; %Shear Modulus
V12 = 0.35;  %Poisson's ration
V23 = 0.35;
disp(eval(C_sol))
disp(eval(inv(C_sol)))

disp(eval(C_sol))