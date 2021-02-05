disp('start symbolic deriviation...')

codegen = ekf_codegen;

%loop period time
syms dt

%process noise convariance matrix
syms Q_i00 Q_i11 Q_i22 Q_i33 Q_i44 Q_i55
codegen = codegen.preload_mat_symbol('Q_i', 6, 6);

%measurement noise covariance matrix
syms V_accel00 V_accel11 V_accel22
codegen = codegen.preload_mat_symbol('V_accel', 3, 3);

syms V_mag00 V_mag11 V_mag22
codegen = codegen.preload_mat_symbol('V_mag', 3, 3);

syms V_gps00 V_gps11 V_gps22 V_gps33
codegen = codegen.preload_mat_symbol('V_gps', 4, 4);

syms V_baro00 V_baro11
codegen = codegen.preload_mat_symbol('V_baro', 2, 2);

%process covariance matrix
syms P_prior00 P_prior01 P_prior02 P_prior03 P_prior04 P_prior05 P_prior06 P_prior07 P_prior08 ...
     P_prior10 P_prior11 P_prior12 P_prior13 P_prior14 P_prior15 P_prior16 P_prior17 P_prior18 ...
     P_prior20 P_prior21 P_prior22 P_prior23 P_prior24 P_prior25 P_prior26 P_prior27 P_prior28 ...
     P_prior30 P_prior31 P_prior32 P_prior33 P_prior34 P_prior35 P_prior36 P_prior37 P_prior38 ...
     P_prior40 P_prior41 P_prior42 P_prior43 P_prior44 P_prior45 P_prior46 P_prior47 P_prior48 ...
     P_prior50 P_prior51 P_prior52 P_prior53 P_prior54 P_prior55 P_prior56 P_prior57 P_prior58 ...
     P_prior60 P_prior61 P_prior62 P_prior63 P_prior64 P_prior65 P_prior66 P_prior67 P_prior68 ...
     P_prior70 P_prior71 P_prior72 P_prior73 P_prior74 P_prior75 P_prior76 P_prior77 P_prior78 ...
     P_prior80 P_prior81 P_prior82 P_prior83 P_prior84 P_prior85 P_prior86 P_prior87 P_prior88
codegen = codegen.preload_mat_symbol('P_prior', 9, 9);

syms P_post00 P_post01 P_post02 P_post03 P_post04 P_post05 P_post06 P_post07 P_post08 ...
     P_post10 P_post11 P_post12 P_post13 P_post14 P_post15 P_post16 P_post17 P_post18 ...
     P_post20 P_post21 P_post22 P_post23 P_post24 P_post25 P_post26 P_post27 P_post28 ...
     P_post30 P_post31 P_post32 P_post33 P_post34 P_post35 P_post36 P_post37 P_post38 ...
     P_post40 P_post41 P_post42 P_post43 P_post44 P_post45 P_post46 P_post47 P_post48 ...
     P_post50 P_post51 P_post52 P_post53 P_post54 P_post55 P_post56 P_post57 P_post58 ...
     P_post60 P_post61 P_post62 P_post63 P_post64 P_post65 P_post66 P_post67 P_post68 ...
     P_post70 P_post71 P_post72 P_post73 P_post74 P_post75 P_post76 P_post77 P_post78 ...
     P_post80 P_post81 P_post82 P_post83 P_post84 P_post85 P_post86 P_post87 P_post88
codegen = codegen.preload_mat_symbol('P_post', 9, 9);

%gamma = sqrt(mx^2 + my^2)
syms gamma

%state vector = position + velocity + quaternion
syms px py pz vx vy vz q0 q1 q2 q3

%magnetometer's reading
syms mx my mz

%accelerometer's reading
syms gx gy gz %gravity's component only

%gps' reading
syms px_gps py_gps vx_gps vy_gps

%accelerometer without bias
syms accel_b_x accel_b_y accel_b_z

%gyroscope without bias
syms gyro_b_x gyro_b_y gyro_b_z

%rotation matrix
syms R00 R01 R02 ...
     R10 R11 R12 ...
     R20 R21 R22
codegen = codegen.preload_mat_symbol('R', 3, 3);

%-R{a_m - a_b}_x * dt
syms R_am_ab_dt00  R_am_ab_dt01  R_am_ab_dt02 ...
     R_am_ab_dt10  R_am_ab_dt11  R_am_ab_dt12 ...
     R_am_ab_dt20  R_am_ab_dt21  R_am_ab_dt22
codegen = codegen.preload_mat_symbol('R_am_ab_dt', 3, 3);

%Rt{(w_m - w_b) * dt}
syms Rt_wm_wb_dt00 Rt_wm_wb_dt01 Rt_wm_wb_dt02 ...
     Rt_wm_wb_dt10 Rt_wm_wb_dt11 Rt_wm_wb_dt12 ...
     Rt_wm_wb_dt20 Rt_wm_wb_dt21 Rt_wm_wb_dt22
codegen = codegen.preload_mat_symbol('Rt_wm_wb_dt', 3, 3);

%kalman gain of accelerometer correction
syms K_accel00 K_accel01 K_accel02 ...
     K_accel10 K_accel11 K_accel12 ...
     K_accel20 K_accel21 K_accel22 ...
     K_accel30 K_accel31 K_accel32 ...
     K_accel40 K_accel41 K_accel42 ...
     K_accel50 K_accel51 K_accel52 ...
     K_accel60 K_accel61 K_accel62 ...
     K_accel70 K_accel71 K_accel72 ...
     K_accel80 K_accel81 K_accel82
codegen = codegen.preload_mat_symbol('K_accel', 9, 3);

%kalman gain of magnetometer correction
syms K_mag00 K_mag01 K_mag02 ...
     K_mag10 K_mag11 K_mag12 ...
     K_mag20 K_mag21 K_mag22 ...
     K_mag30 K_mag31 K_mag32 ...
     K_mag40 K_mag41 K_mag42 ...
     K_mag50 K_mag51 K_mag52 ...
     K_mag60 K_mag61 K_mag62 ...
     K_mag70 K_mag71 K_mag72 ...
     K_mag80 K_mag81 K_mag82
codegen = codegen.preload_mat_symbol('K_mag', 9, 3);

%kalman gain of gps correction
syms K_gps00 K_gps01 K_gps02 Kgps03 ...
     K_gps10 K_gps11 K_gps12 Kgps13 ...
     K_gps20 K_gps21 K_gps22 Kgps23 ...
     K_gps30 K_gps31 K_gps32 Kgps33 ...
     K_gps40 K_gps41 K_gps42 Kgps43 ...
     K_gps50 K_gps51 K_gps52 Kgps53 ...
     K_gps60 K_gps61 K_gps62 Kgps63 ...
     K_gps70 K_gps71 K_gps72 Kgps73 ...
     K_gps80 K_gps81 K_gps82 Kgps83
codegen = codegen.preload_mat_symbol('K_gps', 9, 4);

%============%
% prediction %
%============%

R = [R00 R01 R02;
     R10 R11 R12;
     R20 R21 R22];

skew_a_b = [         0 -accel_b_z +accel_b_y;
            +accel_b_z          0 -accel_b_x;
            -accel_b_y +accel_b_x          0];

R_am_ab_dt = -R * skew_a_b .* dt;

skew_gyro_b = [        0 -gyro_b_z +gyro_b_y;
               +gyro_b_z         0 -gyro_b_x;
               -gyro_b_y +gyro_b_x         0];

Rt_wm_wb_dt = eye(3) - skew_gyro_b .* dt;

%simplfy P_post, assume partial linear independent
Q_i = [[Q_i00     0     0     0     0     0];
       [    0 Q_i11     0     0     0     0];
       [    0     0 Q_i22     0     0     0];
       [    0     0     0 Q_i33     0     0];
       [    0     0     0     0 Q_i44     0];
       [    0     0     0     0     0 Q_i55]];

P_post = ...
    [[P_post00 P_post01 P_post02 P_post03 P_post04 P_post05 P_post06 P_post07 P_post08];
     [P_post10 P_post11 P_post12 P_post13 P_post14 P_post15 P_post16 P_post17 P_post18];
     [P_post20 P_post21 P_post22 P_post23 P_post24 P_post25 P_post26 P_post27 P_post28];
     [P_post30 P_post31 P_post32 P_post33 P_post34 P_post35 P_post36 P_post37 P_post38];
     [P_post40 P_post41 P_post42 P_post43 P_post44 P_post45 P_post46 P_post47 P_post48];
     [P_post50 P_post51 P_post52 P_post53 P_post54 P_post55 P_post56 P_post57 P_post58];
     [P_post60 P_post61 P_post62 P_post63 P_post64 P_post65 P_post66 P_post67 P_post68];
     [P_post70 P_post71 P_post72 P_post73 P_post74 P_post75 P_post76 P_post77 P_post78];
     [P_post80 P_post81 P_post82 P_post83 P_post84 P_post85 P_post86 P_post87 P_post88]];

F_x = [[1 0 0 dt  0  0             0             0             0];
       [0 1 0  0 dt  0             0             0             0];
       [0 0 1  0  0 dt             0             0             0];
       [0 0 0  1  0  0  R_am_ab_dt00  R_am_ab_dt01  R_am_ab_dt02];
       [0 0 0  0  1  0  R_am_ab_dt10  R_am_ab_dt11  R_am_ab_dt12];
       [0 0 0  0  0  1  R_am_ab_dt20  R_am_ab_dt21  R_am_ab_dt22];
       [0 0 0  0  0  0 Rt_wm_wb_dt00 Rt_wm_wb_dt01 Rt_wm_wb_dt02];
       [0 0 0  0  0  0 Rt_wm_wb_dt10 Rt_wm_wb_dt11 Rt_wm_wb_dt12];
       [0 0 0  0  0  0 Rt_wm_wb_dt20 Rt_wm_wb_dt21 Rt_wm_wb_dt22]];

F_i = [[0 0 0 0 0 0];
       [0 0 0 0 0 0];
       [0 0 0 0 0 0];
       [1 0 0 0 0 0];
       [0 1 0 0 0 0];
       [0 0 1 0 0 0];
       [0 0 0 1 0 0];
       [0 0 0 0 1 0];
       [0 0 0 0 0 1]];

P_prior = (F_x * P_post * F_x.') + (F_i * Q_i * F_i.');

%==========================%
% accelerometer correction %
%==========================%

%observation matrix
X_delta_x = [[1 0 0 0 0 0     0     0     0];
             [0 1 0 0 0 0     0     0     0];
             [0 0 1 0 0 0     0     0     0];
             [0 0 0 1 0 0     0     0     0];
             [0 0 0 0 1 0     0     0     0];
             [0 0 0 0 0 1     0     0     0];
             [0 0 0 0 0 0 -q1/2 -q2/2 -q3/2];
             [0 0 0 0 0 0  q0/2 -q3/2  q2/2];
             [0 0 0 0 0 0  q3/2  q0/2 -q1/2];
             [0 0 0 0 0 0 -q2/2  q1/2  q0/2];];

H_x_accel = [[0 0 0 0 0 0  2*q2  2*q3  2*q0 2*q1];
             [0 0 0 0 0 0 -2*q1 -2*q0  2*q3 2*q2];
             [0 0 0 0 0 0  2*q0 -2*q1 -2*q2 2*q3]];

H_accel = H_x_accel * X_delta_x;

%measurement
y_accel = [[gx];
           [gy];
           [gz]];

%observation vector
h_accel = [[2*(q0*q2 + q1*q3)];
           [2*(q2*q3 - q0*q1)];
           [q0*q0 - q1*q1 - q2*q2 + q3*q3]];

resid_accel = y_accel - h_accel;

%measurement noise matrix
V_accel = [[V_accel00         0         0];
           [        0 V_accel11         0];
           [        0         0 V_accel22]];

%kalman gain
Ht_accel = H_accel.';
PHt_accel = P_prior * Ht_accel;
HPHt_V_accel = (H_accel * P_prior * Ht_accel) + V_accel;

K_accel = [[K_accel00 K_accel01 K_accel02];
           [K_accel10 K_accel11 K_accel12];
           [K_accel20 K_accel21 K_accel22];
           [K_accel30 K_accel31 K_accel32];
           [K_accel40 K_accel41 K_accel42];
           [K_accel50 K_accel51 K_accel52];
           [K_accel60 K_accel61 K_accel62];
           [K_accel70 K_accel71 K_accel72];
           [K_accel80 K_accel81 K_accel82]];

%error state innovation
delta_x_accel = K_accel * resid_accel;

%a posteriori covariance matrix update
P_prior_accel = ...
    [[P_prior00 P_prior01 P_prior02 P_prior03 P_prior04 P_prior05 P_prior06 P_prior07 P_prior08];
     [P_prior10 P_prior11 P_prior12 P_prior13 P_prior14 P_prior15 P_prior16 P_prior17 P_prior18];
     [P_prior20 P_prior21 P_prior22 P_prior23 P_prior24 P_prior25 P_prior26 P_prior27 P_prior28];
     [P_prior30 P_prior31 P_prior32 P_prior33 P_prior34 P_prior35 P_prior36 P_prior37 P_prior38];
     [P_prior40 P_prior41 P_prior42 P_prior43 P_prior44 P_prior45 P_prior46 P_prior47 P_prior48];
     [P_prior50 P_prior51 P_prior52 P_prior53 P_prior54 P_prior55 P_prior56 P_prior57 P_prior58];
     [P_prior60 P_prior61 P_prior62 P_prior63 P_prior64 P_prior65 P_prior66 P_prior67 P_prior68];
     [P_prior70 P_prior71 P_prior72 P_prior73 P_prior74 P_prior75 P_prior76 P_prior77 P_prior78];
     [P_prior80 P_prior81 P_prior82 P_prior83 P_prior84 P_prior85 P_prior86 P_prior87 P_prior88]];

P_post_accel = (eye(9) - K_accel*H_accel) * P_prior_accel;

%=========================%
% magnetometer correction %
%=========================%

%observation matrix
H_x_mag = [[0 0 0 0 0 0 2*(+gamma*q0 + mz*q2) 2*(+gamma*q1 + mz*q3) 2*(-gamma*q2 + mz*q0) ...
                                                                    2*(-gamma*q3 + mz*q1)];
           [0 0 0 0 0 0 2*(+gamma*q3 - mz*q1) 2*(+gamma*q2 - mz*q0) 2*(+gamma*q1 + mz*q3) ...
                                                                    2*(+gamma*q0 + mz*q2)];
           [0 0 0 0 0 0 2*(-gamma*q2 + mz*q0) 2*(+gamma*q3 - mz*q1) 2*(-gamma*q0 - mz*q2) ...
                                                                    2*(+gamma*q1 + mz*q3)]];
H_mag = H_x_mag * X_delta_x;

%measurement
y_mag = [[mx];
         [my];
         [mz]];

%observation vector
h_mag = [[gamma*(q0*q0 + q1*q1 - q2*q2 - q3*q3) + 2*mz*(q0*q2 + q1*q3)];
         [2*(gamma*(q1*q2 + q0*q3) + mz*(q2*q3 - q0*q1))];
         [2*gamma*(q1*q3 - q0*q2) + mz*(q0*q0 - q1*q1 - q2*q2 + q3*q3)]];

resid_mag = y_mag - h_mag;

%measurement noise matrix
V_mag = [[V_mag00       0       0];
         [      0 V_mag11       0];
         [      0       0 V_mag22]];

%kalman gain
Ht_mag = H_mag.';
PHt_mag = P_prior * Ht_mag;
HPHt_V_mag = (H_mag * P_prior * Ht_mag) + V_mag;

K_mag = [[K_mag00 K_mag01 K_mag02];
         [K_mag10 K_mag11 K_mag12];
         [K_mag20 K_mag21 K_mag22];
         [K_mag30 K_mag31 K_mag32];
         [K_mag40 K_mag41 K_mag42];
         [K_mag50 K_mag51 K_mag52];
         [K_mag60 K_mag61 K_mag62];
         [K_mag70 K_mag71 K_mag72];
         [K_mag80 K_mag81 K_mag82]];

%error state innovation
delta_x_mag = K_mag * resid_mag;

%a posteriori covariance matrix update
P_prior_mag = ...
    [[P_prior00 P_prior01 P_prior02 P_prior03 P_prior04 P_prior05 P_prior06 P_prior07 P_prior08];
     [P_prior10 P_prior11 P_prior12 P_prior13 P_prior14 P_prior15 P_prior16 P_prior17 P_prior18];
     [P_prior20 P_prior21 P_prior22 P_prior23 P_prior24 P_prior25 P_prior26 P_prior27 P_prior28];
     [P_prior30 P_prior31 P_prior32 P_prior33 P_prior34 P_prior35 P_prior36 P_prior37 P_prior38];
     [P_prior40 P_prior41 P_prior42 P_prior43 P_prior44 P_prior45 P_prior46 P_prior47 P_prior48];
     [P_prior50 P_prior51 P_prior52 P_prior53 P_prior54 P_prior55 P_prior56 P_prior57 P_prior58];
     [P_prior60 P_prior61 P_prior62 P_prior63 P_prior64 P_prior65 P_prior66 P_prior67 P_prior68];
     [P_prior70 P_prior71 P_prior72 P_prior73 P_prior74 P_prior75 P_prior76 P_prior77 P_prior78];
     [P_prior80 P_prior81 P_prior82 P_prior83 P_prior84 P_prior85 P_prior86 P_prior87 P_prior88]];

P_post_mag = (eye(9) - K_mag*H_mag) * P_prior_mag;

%================%
% gps correction %
%================%

%observation matrix
H_x_gps = [[1 0 0 0 0 0 0 0 0 0];
           [0 1 0 0 0 0 0 0 0 0];
           [0 0 0 1 0 0 0 0 0 0];
	   [0 0 0 0 1 0 0 0 0 0]];
H_gps = H_x_gps * X_delta_x;

%measurement (from gps receiver)
y_gps = [[px_gps];
         [py_gps];
         [vx_gps];
	 [vy_gps]];

%observation vector (from state vector)
h_gps = [[px];
         [py];
         [vx];
	 [vy]];

resid_gps = y_gps - h_gps;

%measurement noise matrix
V_gps = [[V_gps00       0       0       0];
         [      0 V_gps11       0       0];
         [      0       0 V_gps22       0];
	 [      0       0       0 V_gps33]];

%kalman gain
Ht_gps = H_gps.';
PHt_gps = P_prior * Ht_gps;
HPHt_V_gps = (H_gps * P_prior * Ht_gps) + V_gps;

K_gps = [[K_gps00 K_gps01 K_gps02 Kgps03];
         [K_gps10 K_gps11 K_gps12 Kgps13];
         [K_gps20 K_gps21 K_gps22 Kgps23];
         [K_gps30 K_gps31 K_gps32 Kgps33];
         [K_gps40 K_gps41 K_gps42 Kgps43];
         [K_gps50 K_gps51 K_gps52 Kgps53];
         [K_gps60 K_gps61 K_gps62 Kgps63];
         [K_gps70 K_gps71 K_gps72 Kgps73];
         [K_gps80 K_gps81 K_gps82 Kgps83]];

%error state innovation
delta_x_gps = K_gps * resid_gps;

%a posteriori covariance matrix update
P_prior_gps = ...
    [[P_prior00 P_prior01 P_prior02 P_prior03 P_prior04 P_prior05 P_prior06 P_prior07 P_prior08];
     [P_prior10 P_prior11 P_prior12 P_prior13 P_prior14 P_prior15 P_prior16 P_prior17 P_prior18];
     [P_prior20 P_prior21 P_prior22 P_prior23 P_prior24 P_prior25 P_prior26 P_prior27 P_prior28];
     [P_prior30 P_prior31 P_prior32 P_prior33 P_prior34 P_prior35 P_prior36 P_prior37 P_prior38];
     [P_prior40 P_prior41 P_prior42 P_prior43 P_prior44 P_prior45 P_prior46 P_prior47 P_prior48];
     [P_prior50 P_prior51 P_prior52 P_prior53 P_prior54 P_prior55 P_prior56 P_prior57 P_prior58];
     [P_prior60 P_prior61 P_prior62 P_prior63 P_prior64 P_prior65 P_prior66 P_prior67 P_prior68];
     [P_prior70 P_prior71 P_prior72 P_prior73 P_prior74 P_prior75 P_prior76 P_prior77 P_prior78];
     [P_prior80 P_prior81 P_prior82 P_prior83 P_prior84 P_prior85 P_prior86 P_prior87 P_prior88]];

P_post_gps = (eye(9) - K_gps*H_gps) * P_prior_gps;

%========================%
% save derivation result %
%========================%

tic();

disp('start code generation...')

%============%
% prediction %
%============%
codegen = codegen.open_file('predict.c');

codegen.add_c_comment('/* calculate a priori process covariance matrix */');
codegen.generate_c_code('R_am_ab_dt', R_am_ab_dt, 'is_symmetry=0');
codegen.generate_c_code('Rt_wm_wb_dt', Rt_wm_wb_dt, 'is_symmetry=0');
codegen.generate_c_code('P_prior', P_prior, 'is_symmetry=1')

codegen.close_file();

%==========================%
% accelerometer correction %
%==========================%
codegen = codegen.open_file('accelerometer_correct.c');

codegen.add_c_comment('/* calculate kalman gain subterm P * transpose(H) */');
codegen.generate_c_code('PHt_accel', PHt_accel, 'is_symmetry=0')
%codegen.generate_c_code('HPHt_V_accel', HPHt_V_accel)

codegen.add_c_comment('/* calculate error state residual */');
codegen.generate_c_code('delta_x_accel', delta_x_accel, 'is_symmetry=0')

codegen.add_c_comment('/* calculate a posteriori process covariance matrix */');
codegen.generate_c_code('P_post_accel', P_post_accel, 'is_symmetry=1')

codegen.close_file();

%=========================%
% magnetometer correction %
%=========================%
codegen = codegen.open_file('magnetometer_correct.c');

codegen.add_c_comment('/* calculate kalman gain subterm P * transpose(H) */');
codegen.generate_c_code('PHt_mag', PHt_gps, 'is_symmetry=0')
%codegen.generate_c_code('HPHt_V_mag', HPHt_V_mag)

codegen.add_c_comment('/* calculate error state residual */');
codegen.generate_c_code('delta_x_mag', delta_x_gps, 'is_symmetry=0')

codegen.add_c_comment('/* calculate a posteriori process covariance matrix */');
codegen.generate_c_code('P_post_mag', P_post_gps, 'is_symmetry=1')

codegen.close_file();

%================%
% gps correction %
%================%
codegen = codegen.open_file('gps.c');

codegen.add_c_comment('/* calculate kalman gain subterm P * transpose(H) */');
codegen.generate_c_code('PHt_gps', PHt_mag, 'is_symmetry=0')
%codegen.generate_c_code('HPHt_V_mag', HPHt_V_mag)

codegen.add_c_comment('/* calculate error state residual */');
codegen.generate_c_code('delta_x_gps', delta_x_mag, 'is_symmetry=0')

codegen.add_c_comment('/* calculate a posteriori process covariance matrix */');
codegen.generate_c_code('P_post_gps', P_post_mag, 'is_symmetry=1')

codegen.close_file();

toc();
