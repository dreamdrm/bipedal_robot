clc; clear; close all;

quy_dao_3;

q1_sol = double(danh_sach_nghiem.q1(1));
q2_sol = double(danh_sach_nghiem.q2(1));
q3_sol = double(danh_sach_nghiem.q3(1));

syms dq1 dq2 dq3 real;
q = q_syms; 
dq = [dq1; dq2; dq3];

syms d1 a3 g real
syms m1 m2 m3 real
syms xc1 yc1 zc1 xc2 yc2 zc2 xc3 yc3 zc3 real
syms Ixx1 Iyy1 Izz1 Ixx2 Iyy2 Izz2 Ixx3 Iyy3 Izz3 real

gia_tri_so = [
    
    g, 9.81;
    
    % khoi luong
    m1, 2.0; 
    m2, 2.2; 
    m3, 0.75;
    
    % toa do khoi tam
    xc1, 0.0;        yc1, 0.12088;    zc1, 0.0;
    xc2, 0.0;        yc2, 0.22369;    zc2, -0.00209;
    xc3, -0.10125;   yc3, 0.0;        zc3, 0.0;
    
    % tensor quán tính chính (kg.m^2)
    Ixx1, 0.013876;  Iyy1, 0.002864;  Izz1, 0.013872;
    Ixx2, 0.042514;  Iyy2, 0.001610;  Izz2, 0.042724;
    Ixx3, 0.000458;  Iyy3, 0.003788;  Izz3, 0.003831;
];

% ma tran quán tính chính
I1 = diag([Ixx1, Iyy1, Izz1]);
I2 = diag([Ixx2, Iyy2, Izz2]);
I3 = diag([Ixx3, Iyy3, Izz3]);

% Vector toa do thuan nhat
Pc1_1 = [xc1; yc1; zc1; 1];
Pc2_2 = [xc2; yc2; zc2; 1];
Pc3_3 = [xc3; yc3; zc3; 1];


H01 = H_all{1};
H02 = H01 * H_all{2};
H03 = H02 * H_all{3};

% ma tran quay
R1 = H01(1:3, 1:3); R2 = H02(1:3, 1:3); R3 = H03(1:3, 1:3);

% Vector toa do decartes
Pc1_0 = H01 * Pc1_1; r_c1_0 = Pc1_0(1:3);
Pc2_0 = H02 * Pc2_2; r_c2_0 = Pc2_0(1:3);
Pc3_0 = H03 * Pc3_3; r_c3_0 = Pc3_0(1:3);

% ma tra jacobi tinh tien
J_c1 = jacobian(r_c1_0, q);
J_c2 = jacobian(r_c2_0, q);
J_c3 = jacobian(r_c3_0, q);

dR1 = sym(zeros(3,3)); dR2 = sym(zeros(3,3)); dR3 = sym(zeros(3,3));
for i = 1:3
    dR1 = dR1 + diff(R1, q(i)) * dq(i);
    dR2 = dR2 + diff(R2, q(i)) * dq(i);
    dR3 = dR3 + diff(R3, q(i)) * dq(i);
end
S1 = dR1 * R1.'; omega1 = [S1(3,2); S1(1,3); S1(2,1)];
S2 = dR2 * R2.'; omega2 = [S2(3,2); S2(1,3); S2(2,1)];
S3 = dR3 * R3.'; omega3 = [S3(3,2); S3(1,3); S3(2,1)];

% ma tran jacobi quay
Jw1 = jacobian(omega1, dq);
Jw2 = jacobian(omega2, dq);
Jw3 = jacobian(omega3, dq);

m = [m1, m2, m3];
I_cell  = {I1, I2, I3};
R_cell  = {R1, R2, R3};
J_c_cell = {J_c1, J_c2, J_c3};
Jw_cell = {Jw1, Jw2, Jw3};
n = 3;
M_sym = sym(zeros(n, n));

% ma tran khoi luong
for i = 1:n
    M_tinh_tien = J_c_cell{i}.'* m(i) * J_c_cell{i};
  
    M_quay = Jw_cell{i}.' * (R_cell{i} * I_cell{i} * R_cell{i}.') * Jw_cell{i};
    
    M_sym = M_sym + M_tinh_tien + M_quay;
end

M_sym = simplify(M_sym);
M_sym = subs(M_sym, [d1, a3], [d1_val, a3_val]);

% Vector trong luc
P_sym = g * (m1*r_c1_0(3) + m2*r_c2_0(3) + m3*r_c3_0(3));
G_sym = jacobian(P_sym, q).';
G_sym = subs(G_sym, [d1, a3], [d1_val, a3_val]);

% Ma tran Coriolis
C_sym = sym(zeros(n, n));
for k = 1:n
    for j = 1:n
        for i = 1:n
            c_ijk = 0.5 * (diff(M_sym(k,j), q(i)) + diff(M_sym(k,i), q(j)) - diff(M_sym(i,j), q(k)));
            C_sym(k,j) = C_sym(k,j) + c_ijk * dq(i);
        end
    end
end

vars = gia_tri_so(:, 1);
vals = gia_tri_so(:, 2);

M_sym = subs(M_sym, vars, vals);
G_sym = subs(G_sym, vars, vals);
C_sym = subs(C_sym, vars, vals);

disp('MA TRAN KHOI LUONG M(q)'); 
disp(vpa(M_sym, 4));

disp('MA TRAN CORIOLIS C(q, dq)'); 
disp(vpa(C_sym, 4));

disp('VECTOR TRONG LUC G(q)'); 
disp(vpa(G_sym, 4));


calc_M = matlabFunction(M_sym, 'Vars', {q});
calc_G = matlabFunction(G_sym, 'Vars', {q});
calc_C = matlabFunction(C_sym, 'Vars', {q, dq});

num_points = length(t);
tau = zeros(3, num_points); 

for k = 1:num_points
    
    q_k   = Q(:, k);
    dq_k  = dQ(:, k);
    ddq_k = ddQ(:, k);
    
    M_k = calc_M(q_k);
    C_k = calc_C(q_k, dq_k);
    G_k = calc_G(q_k);
  
    tau(:, k) = M_k * ddq_k + C_k * dq_k + G_k;
end

figure('Name', 'BIEU DO LUC VA MOMEN', 'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);

% khop 1 (N.m)
subplot(3, 1, 1);
plot(t, tau(1, :), 'r', 'LineWidth', 2);
grid on; ylabel('Momen \tau_1 (N.m)');
title('Momen yeu cau cho khop quay 1');

% khop 2 (N)
subplot(3, 1, 2);
plot(t, tau(2, :), 'g', 'LineWidth', 2);
grid on; ylabel('Luc \tau_2 (N)');
title('Luc day yeu cau cho khop tinh tien 2');

% khop 3 (N.m)
subplot(3, 1, 3);
plot(t, tau(3, :), 'b', 'LineWidth', 2);
grid on; xlabel('Thoi gian (s)'); ylabel('Momen \tau_3 (N.m)');
title('Momen yeu cau cho khop quay 3');

tau_max = [max(abs(tau(1, :))); 
           max(abs(tau(2, :))); 
           max(abs(tau(3, :)))];

fprintf('Thông s? công su?t c?c ??i c?n thi?t (Max Capacity):\n');
fprintf('Kh?p 1 (Quay)       : Tau_1_max = %.4f (N.m)\n', tau_max(1));
fprintf('Kh?p 2 (T?nh ti?n)  : Tau_2_max = %.4f (N)\n', tau_max(2));
fprintf('Kh?p 3 (Quay)       : Tau_3_max = %.4f (N.m)\n\n', tau_max(3));

% 2. Gi? l?p c?p l?c max vào v? trí xu?t phát (t = 0)
q_test  = Q(:, 1);
dq_test = dQ(:, 1); 

% Tính các ma tr?n h? th?ng t?i t? th? xu?t phát
M_test = calc_M(q_test);
C_test = calc_C(q_test, dq_test);
G_test = calc_G(q_test);

% 3. Gi?i ph??ng trình ??ng l?c h?c thu?n: ddq = M \ (tau - C*dq - G)
ddq_thuan = M_test \ (tau_max - C_test * dq_test - G_test);

disp('=> Gia t?c t?c th?i l?n nh?t sinh ra khi c?p l?c c?c ??i t?i v? trí xu?t phát:');
fprintf('Gia t?c kh?p 1 (ddq_1) = %8.4f (rad/s^2)\n', ddq_thuan(1));
fprintf('Gia t?c kh?p 2 (ddq_2) = %8.4f (m/s^2)\n', ddq_thuan(2));
fprintf('Gia t?c kh?p 3 (ddq_3) = %8.4f (rad/s^2)\n', ddq_thuan(3));