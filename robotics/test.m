% =========================================================================
% CH??NG TRÌNH TÍNH ??NG H?C THU?N CÓ KH?P T?NH TI?N
% =========================================================================
clear; clc;

% Khai báo các bi?n ký hi?u 
% q1, q3: bi?n góc quay (theta)
% q2: bi?n t?nh ti?n (d)
syms q1 q2 q3 l1 l3 real 

% B?ng thông s? D-H (C?u trúc: [theta, d, a, alpha])
% L?u ý: C?t 2 c?a Khâu 2 bây gi? ch?a bi?n t?nh ti?n q2. 
% Các h?ng s? (0, pi/2, -pi/2...) d??i ?ây là ví d?, b?n hãy ??i chi?u 
% l?i v?i b?ng D-H b?n l?p t? s? ?? ?? thay s? cho chính xác nhé!
DH_Table = [
    q1,    l1,    0,   pi/2;  % Khâu 1: Kh?p quay (bi?n q1 n?m ? c?t theta)
     0,    q2,    0,  -pi/2;  % Khâu 2: Kh?p T?NH TI?N (bi?n q2 n?m ? c?t d)
    q3,     0,   l3,      0   % Khâu 3: Kh?p quay (bi?n q3 n?m ? c?t theta)
];

num_links = size(DH_Table, 1);
T_0_n = eye(4); 

fprintf('--- ?ANG TÍNH TOÁN CÁC MA TR?N BI?N ??I ---\n\n');

for i = 1:num_links
    theta = DH_Table(i, 1);
    d     = DH_Table(i, 2);
    a     = DH_Table(i, 3);
    alpha = DH_Table(i, 4);
    
    % Công th?c ma tr?n bi?n ??i D-H chu?n
    A = [ cos(theta), -sin(theta)*cos(alpha),  sin(theta)*sin(alpha), a*cos(theta);
          sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
                   0,             sin(alpha),             cos(alpha),           d;
                   0,                      0,                      0,           1 ];
    
    fprintf('Ma tr?n A_%d_%d:\n', i-1, i);
    disp(simplify(A));
    
    T_0_n = T_0_n * A;
end

T_0_n = simplify(T_0_n);

disp('=========================================================================');
disp('MA TR?N BI?N ??I THU?N NH?T T?NG QUÁT (T_0_n):');
disp(T_0_n);

% V? trí ?i?m tác ??ng cu?i
x_E = T_0_n(1, 4);
y_E = T_0_n(2, 4);
z_E = T_0_n(3, 4);

disp('PH??NG TRÌNH V? TRÍ T?A ?? BÀN TAY MÁY:');
fprintf('Px = %s\n\n', char(x_E));
fprintf('Py = %s\n\n', char(y_E));
fprintf('Pz = %s\n\n', char(z_E));
