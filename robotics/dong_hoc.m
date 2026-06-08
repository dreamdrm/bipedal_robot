function [T_0_n, P_E, H_all, R_0_i_all, danh_sach_nghiem, DH_Table, d1_val, a3_val, q_syms] = dong_hoc()

    syms q1 d1 q2 q3 a3 real; 
    q_syms = [q1; q2; q3];
    
    DH_Table = [
        q1-pi/2, d1, 0, -pi/2;
        pi/2,    q2, 0, -pi/2;
        -q3,     0,  a3,  0
        ];

    disp('BANG DH:');
    disp(DH_Table);
    num_links = size(DH_Table, 1);

    T_0_n = eye(4); 
    H_all = cell(1, num_links);
    R_0_i_all = cell(1, num_links);
    
    for i = 1:num_links
        theta = DH_Table(i, 1);
        d     = DH_Table(i, 2);
        a     = DH_Table(i, 3);
        alpha = DH_Table(i, 4);


        H = [ cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha) a*cos(theta);
              sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta);
                       0             sin(alpha)             cos(alpha)           d;
                       0                      0                      0           1 
            ];
        
        H_all{i} = simplify(H);
        fprintf('H_%d_%d:\n', i-1, i);
        disp(H_all{i});

        T_0_n = T_0_n * H;
        
        R_0_i_all{i} = simplify(T_0_n(1:3, 1:3));
        fprintf('MA TRAN QUAY R_0_%d:\n', i);
        disp(R_0_i_all{i});
    end

    T_0_n = simplify(T_0_n);
    

    fprintf('MA TRAN TRANG THAI KHAU CUOI (T_0_%d):\n', num_links);
    disp(T_0_n);

    x_E = T_0_n(1, 4);
    y_E = T_0_n(2, 4);
    z_E = T_0_n(3, 4);
    P_E = [x_E; y_E; z_E];

    disp('TOA DO DIEM CUOI:');
    pretty(P_E);

    Xd = 0;
    Yd = 0.65;
    Zd = 0.25;

    fprintf('x = %.4f\n', Xd);
    fprintf('y = %.4f\n', Yd);
    fprintf('z = %.4f\n\n', Zd);

    d1_val = 0.25;
    a3_val = 0.25;

    Px_eq = subs(x_E, [d1, a3], [d1_val, a3_val]);
    Py_eq = subs(y_E, [d1, a3], [d1_val, a3_val]);
    Pz_eq = subs(z_E, [d1, a3], [d1_val, a3_val]);

    eq1 = Px_eq == Xd;
    eq2 = Py_eq == Yd;
    eq3 = Pz_eq == Zd;


    gioi_han_may = [-pi pi; 0.15 0.4; 0 pi];

    nghiem = vpasolve([eq1, eq2, eq3], [q1, q2, q3], gioi_han_may);
    danh_sach_nghiem = nghiem;
    if isempty(nghiem.q1)
        disp('KHONG TIM THAY NGHIEM!');
    else
        disp('BO NGHIEM TIM DUOC:');

        for k = 1:length(nghiem.q1)

            fprintf('  Goc khop q1      = %.4f (rad) = %.2f°\n', double(nghiem.q1(k)), double(nghiem.q1(k))*180/pi);
            fprintf('  Do tinh tien q2  = %.4f (m)\n', double(nghiem.q2(k)));
            fprintf('  Goc khop q3      = %.4f (rad) = %.2f°\n\n', double(nghiem.q3(k)), double(nghiem.q3(k))*180/pi);
        end
    end
end    
