
[~, ~, H_all, ~, danh_sach_nghiem, ~, d1_val, a3_val, q_syms] = dong_hoc();

if isempty(danh_sach_nghiem.q1)
    error('NGOAI PHAM VI CHUYEN DONG CUA MAY');
end

q_start = [0; 0.15; 0];

q_end = [double(danh_sach_nghiem.q1(1));
         double(danh_sach_nghiem.q2(1));
         double(danh_sach_nghiem.q3(1))];

fprintf('\nXUAT PHAT (Point A): [%.4f, %.4f, %.4f]^T\n', q_start);
fprintf('DIEM DEN (Point B):  [%.4f, %.4f, %.4f]^T\n\n', q_end);

t0 = 0;             
tf = 5;            
t = linspace(t0, tf, 200); 

v_start = [0; 0; 0];
v_end   = [0; 0; 0];

Q   = zeros(3, length(t)); 
dQ  = zeros(3, length(t)); 
ddQ = zeros(3, length(t)); 

for i = 1:3
 
    % q(t) = a0 + a1*t + a2*t^2 + a3*t^3 
    M = [1,  t0,  t0^2,   t0^3;
         0,  1,   2*t0,   3*t0^2;      
         1,  tf,  tf^2,   tf^3;
         0,  1,   2*tf,   3*tf^2;
         ];
    
    b = [q_start(i); v_start(i); q_end(i); v_end(i)];
    
    a = M \ b; 
    
    Q(i,:)   = a(1) + a(2)*t + a(3)*t.^2 + a(4)*t.^3;
    dQ(i,:)  = a(2) + 2*a(3)*t + 3*a(4)*t.^2;
    ddQ(i,:) = 2*a(3) + 6*a(4)*t;
end

figure('Name', 'Quy dao khong gian khop (Vi trí - Van toc - Gia toc)', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
joint_names = {'Khop 1 (q_1) [rad]', 'Khop 2 (q_2) [m]', 'Khop 3 (q_3) [rad]'};

for i = 1:3
    
    % Vi tri
    subplot(3, 3, (i-1)*3 + 1);
    plot(t, Q(i,:), 'b', 'LineWidth', 2);
    grid on; title(['Vi tri ', joint_names{i}]);
    xlabel('Thoi gian t (s)'); ylabel('Vi tri');

    % van toc
    subplot(3, 3, (i-1)*3 + 2);
    plot(t, dQ(i,:), 'g', 'LineWidth', 2);
    grid on; title(['Van toc ', joint_names{i}, '/s']);
    xlabel('Thoi gian t (s)'); ylabel('Van toc');
    
    % gia toc
    subplot(3, 3, (i-1)*3 + 3);
    plot(t, ddQ(i,:), 'r', 'LineWidth', 2);
    grid on; title(['Gia toc ', joint_names{i}, '/s^2']);
    xlabel('Thoi gian t (s)'); ylabel('Gia toc');
end


X_eff = zeros(1, length(t));
Y_eff = zeros(1, length(t));
Z_eff = zeros(1, length(t));

for k = 1:length(t)
    q1_k = Q(1,k);
    q2_k = Q(2,k);
    q3_k = Q(3,k);
    
    H1 = func_DH(q1_k - pi/2, d1_val,   0,  -pi/2);
    H2 = func_DH(pi/2,        q2_k, 0,  -pi/2);
    H3 = func_DH(-q3_k,       0,    a3_val, 0);
    
    T_end = H1 * H2 * H3; 
    
    X_eff(k) = T_end(1,4);
    Y_eff(k) = T_end(2,4);
    Z_eff(k) = T_end(3,4);
end

figure('Name', 'QUY DAO CHUYEN DONG TRONG KHONG GIAN THAO TAC', 'NumberTitle', 'off');
plot3(X_eff, Y_eff, Z_eff, 'k', 'LineWidth', 2.5); hold on; grid on;
plot3(X_eff(1), Y_eff(1), Z_eff(1), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot3(X_eff(end), Y_eff(end), Z_eff(end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

xlabel('Truc X (m)'); ylabel('Truc Y (m)'); zlabel('Truc Z (m)');
title('QUY DAO CHUYEN DONG DIEM CUOI');
legend('Duong di quy dao', 'Diem xuat phat (A)', 'Diem dich (B)');
view(3); 

function H = func_DH(theta, d, a, alpha)
    H = [ cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha) a*cos(theta);
          sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta);
                   0             sin(alpha)             cos(alpha)           d;
                   0                      0                      0           1 ];
end