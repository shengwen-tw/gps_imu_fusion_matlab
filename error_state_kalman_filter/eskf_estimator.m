classdef eskf_estimator
    %assumption: NED coordinate system
    
    properties
        g_constant = 9.8;
        
        home_longitude = 0;
        home_latitude = 0;
        home_ecef_x = 0;
        home_ecef_y = 0;
        home_ecef_z = 0;
        
        EQUATORIAL_RADIUS = 6378137 %[m], earth semi-major length (AE)
        POLAR_RADIUS = 6356752      %[m], earth semi-minor length (AP)
        AP_SQ_DIV_AE_SQ = 0.99331   %(AP^2)/(AE^2)
        ECCENTRICITY = 0.081820     %e^2 = 1 - (AP^2)/(AE^2)
        
        gravity = [0; 0; 9.8];
        
        %nomnial state
        x_nominal = [0;  %px
                     0;  %py
                     0;  %pz
                     0;  %vx
                     0;  %vy
                     0;  %vz
                     1;  %q0
                     0;  %q2
                     0;  %q2
                     0;  %q3
                     0;  %a_b_x
                     0;  %a_b_y
                     0;  %a_b_z
                     0;  %w_b_x
                     0;  %w_b_y
                     0]; %w_b_z
                      
        %error state
        delta_x = [0;  %px
                   0;  %py
                   0;  %pz
                   0;  %vx
                   0;  %vy
                   0;  %vz
                   0;  %theta_x
                   0;  %theta_y
                   0;  %theta_z
                   0;  %a_b_x
                   0;  %a_b_y
                   0;  %a_b_z
                   0;  %w_b_x
                   0;  %w_b_y
                   0]; %w_b_z
                      
        %attitude direction cosine matrix
        R = eye(3);  %body-fixed frame to inertial frame
        Rt = eye(3); %inertial frame to body-fixed frame
        
        %covariance matrix of process white noise
        Q_i = [1e-5 0 0 0 0 0 0 0 0 0 0 0;   %noise of ax
               0 1e-5 0 0 0 0 0 0 0 0 0 0;   %noise of ay
               0 0 1e-5 0 0 0 0 0 0 0 0 0;   %noise of az
               0 0 0 1e-5 0 0 0 0 0 0 0 0;   %noise of wx
               0 0 0 0 1e-5 0 0 0 0 0 0 0;   %noise of wy
               0 0 0 0 0 1e-5 0 0 0 0 0 0;   %noise of wz
               0 0 0 0 0 0 1e-9 0 0 0 0 0;   %noise of a_b_x
               0 0 0 0 0 0 0 1e-9 0 0 0 0;   %noise of a_b_y
               0 0 0 0 0 0 0 0 1e-9 0 0 0;   %noise of a_b_z  
               0 0 0 0 0 0 0 0 0 1e-10 0 0;  %noise of w_b_x
               0 0 0 0 0 0 0 0 0 0 1e-10 0;  %noise of w_b_y
               0 0 0 0 0 0 0 0 0 0 0 1e-10]; %noise of w_b_z
        
        %process covariance matrix of error state
        P = [5 0 0 0 0 0 0 0 0 0 0 0 0 0 0;     %delta px
             0 5 0 0 0 0 0 0 0 0 0 0 0 0 0;     %delta py
             0 0 5 0 0 0 0 0 0 0 0 0 0 0 0;     %delta pz
             0 0 0 5 0 0 0 0 0 0 0 0 0 0 0;     %delta vx
             0 0 0 0 5 0 0 0 0 0 0 0 0 0 0;     %delta vy
             0 0 0 0 0 5 0 0 0 0 0 0 0 0 0;     %delta vz
             0 0 0 0 0 0 5 0 0 0 0 0 0 0 0;     %delta_x
             0 0 0 0 0 0 0 5 0 0 0 0 0 0 0;     %delta_y
             0 0 0 0 0 0 0 0 5 0 0 0 0 0 0;     %delta_z
             0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;     %delta a_b_x
             0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;     %delta a_b_y
             0 0 0 0 0 0 0 0 0 0 0 1e-2 0 0 0;  %delta a_b_z
             0 0 0 0 0 0 0 0 0 0 0 0 1e-1 0 0;  %delta w_b_x
             0 0 0 0 0 0 0 0 0 0 0 0 0 1e-1 0;  %delta w_b_y
             0 0 0 0 0 0 0 0 0 0 0 0 0 0 1e-1]; %delta w_b_z
        
        %observation covariance matrix of accelerometer
        V_accel = [7e-1 0 0;  %ax
                   0 7e-1 0;  %ay
                   0 0 7e-1]; %az
               
        %observation covariance matrix of accelerometer
        V_mag = [5e-1 0 0;  %mx
                 0 5e-1 0;  %my
                 0 0 5e-1]; %mz
             
        %observation covariance matrix of the gps sensor
        V_gps = [1e-4 0 0 0;     %px
                 0 1e-4 0 0;     %py
                 0 0 5e-4 0;  %vx
                 0 0 0 5e-4]; %vy
             
        %%observation covariance matrix of the height sensor
        V_height = [1e-1 0;  %pz
                    0 1e-1]; %vz   
                
        I_3x3 = eye(3);
        I_4x4 = eye(4);
        I_6x6 = eye(6);
        I_9x9 = eye(9);
        I_12x12 = eye(12);
        I_15x15 = eye(15);
    end
    
    methods
        function q_out = quaternion_mult(obj, q1, q2)
            q_out = [q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) - q1(4)*q2(4);
                     q1(1)*q2(2) + q1(2)*q2(1) + q1(3)*q2(4) - q1(4)*q2(3);
                     q1(1)*q2(3) - q1(2)*q2(4) + q1(3)*q2(1) + q1(4)*q2(2);
                     q1(1)*q2(4) + q1(2)*q2(3) - q1(3)*q2(2) + q1(4)*q2(1)];
        end
        
        function q_out = quaternion_conj(obj, q)
            q_out = [+q(1);
                     -q(2);
                     -q(3);
                     -q(4)];
        end
        
        function q_out = quat_normalize(obj, q)
            div_q_norm = 1 / norm(q);
            q_out = [q(1) * div_q_norm;
                     q(2) * div_q_norm;
                     q(3) * div_q_norm;
                     q(4) * div_q_norm];
        end
        
        function R = quat_to_rotation_matrix(obj, q)
            q1q1 = q(2) * q(2);
            q2q2 = q(3) * q(3);
            q3q3 = q(4) * q(4);
            q1q2 = q(2) * q(3);
            q0q2 = q(1) * q(3);
            q0q3 = q(1) * q(4);
            q1q3 = q(2) * q(4);
            q2q3 = q(3) * q(4);
            q0q1 = q(1) * q(2);

            r11 = 1.0 - 2.0 * (q2q2 + q3q3);
            r12 = 2.0 * (q1q2 - q0q3);
            r13 = 2.0 * (q0q2 + q1q3);

            r21 = 2.0 * (q1q2 + q0q3);
            r22 = 1.0 - 2.0 * (q1q1 + q3q3);
            r23 = 2.0 * (q2q3 - q0q1);

            r31 = 2.0 * (q1q3 - q0q2);
            r32 = 2.0 * (q0q1 + q2q3);
            r33 = 1.0 - 2.0 * (q1q1 + q2q2);
            
            R = [r11 r12 r13;
                 r21 r22 r23;
                 r31 r32 r33];
        end

        function euler_angles = quat_to_euler(obj, q)
            roll = atan2(2.0*(q(1)*q(2) + q(3)*q(4)), 1.0-2.0*(q(2)*q(2) + q(3)*q(3)));
            pitch = asin(2.0*(q(1)*q(3) - q(4)*q(2)));
            yaw = atan2(2.0*(q(1)*q(4) + q(2)*q(3)), 1.0-2.0*(q(3)*q(3) + q(4)*q(4)));
            
            roll = rad2deg(roll);
            pitch = rad2deg(pitch);
            yaw = rad2deg(yaw);
            
            euler_angles = [roll; pitch; yaw];
        end
        
        function skew_matrix = hat_map_3x3(obj, skew_vector)
            skew_matrix(1, 1) = 0;
            skew_matrix(1, 2) = -skew_vector(3);
            skew_matrix(1, 3) = +skew_vector(2);
            skew_matrix(2, 1) = +skew_vector(3);
            skew_matrix(2, 2) = 0;
            skew_matrix(2, 3) = -skew_vector(1);
            skew_matrix(3, 1) = -skew_vector(2);
            skew_matrix(3, 2) = +skew_vector(1);
            skew_matrix(3, 3) = 0;
        end
        
        function quaternion = get_quaternion(obj)
            quaternion = obj.x_nominal(7:10);
        end
        
        function vec_enu = convert_3x1_vector_ned_to_enu(obj, vec_ned)
            vec_enu = [+vec_ned(2);
                       +vec_ned(1);
                       -vec_ned(3)];
        end
        
        function ret_obj = set_home_longitude_latitude(obj, longitude, latitude, height_msl)
            obj.home_longitude = longitude;
            obj.home_latitude = latitude;
            
            sin_lambda = sin(deg2rad(longitude));
            cos_lambda = cos(deg2rad(longitude));
            sin_phi = sin(deg2rad(latitude));
            cos_phi = cos(deg2rad(latitude));
            
            %convert geodatic coordinates to earth center earth fixed frame (ecef)
            N = obj.EQUATORIAL_RADIUS / sqrt(1 - (obj.ECCENTRICITY * sin_phi * sin_phi));
            obj.home_ecef_x = (N + height_msl) * cos_phi * cos_lambda;
            obj.home_ecef_y = (N + height_msl) * cos_phi * sin_lambda;
            obj.home_ecef_z = (obj.AP_SQ_DIV_AE_SQ * N + height_msl) * sin_phi;
            
            ret_obj = obj;
        end
        
        function enu_pos = covert_geographic_to_ned_frame(obj, longitude, latitude, height_msl)
            sin_lambda = sin(deg2rad(longitude));
            cos_lambda = cos(deg2rad(longitude));
            sin_phi = sin(deg2rad(latitude));
            cos_phi = cos(deg2rad(latitude));

            %convert geodatic coordinates to earth center earth fixed frame (ecef)
            N = obj.EQUATORIAL_RADIUS / sqrt(1 - (obj.ECCENTRICITY * sin_phi * sin_phi));
            ecef_now_x = (N + height_msl) * cos_phi * cos_lambda;
            ecef_now_y = (N + height_msl) * cos_phi * sin_lambda;
            ecef_now_z = (obj.AP_SQ_DIV_AE_SQ * N + height_msl) * sin_phi;
            
            %convert position from earth center earth fixed frame to east north up frame
            r11 = -sin_phi * cos_lambda;
            r12 = -sin_lambda;
            r13 = -cos_phi * cos_lambda;
            r21 = -sin_phi * sin_lambda;
            r22 = cos_lambda;
            r23 = -cos_phi * sin_lambda;
            r31 = cos_phi;
            r32 = 0;
            r33 = -sin_phi;
            R = [r11, r12, r13;
                 r21, r22, r23
                 r31, r32, r33];
            
            dx = ecef_now_x - obj.home_ecef_x;
            dy = ecef_now_y - obj.home_ecef_y;
            dz = ecef_now_z - obj.home_ecef_z;
            
            enu_pos = R.' * [dx; dy; dz];
        end
        
        %for initialization only
        function q = convert_gravity_to_quat(obj, a)
            if a(3) >= 0.0
                sqrt_tmp = 1.0 / sqrt(2.0 * (a(3) + 1));
                q(1) = sqrt(0.5 * (a(3) + 1.0));
                q(2) = -a(2) * sqrt_tmp;
                q(3) = +a(1) * sqrt_tmp;
                q(4) = 0.0;
            else
                sqrt_tmp = 1.0 / sqrt(2.0 * (1.0 - a(3)));
                q(1) = -a(2) * sqrt_tmp;
                q(2) = sqrt((1.0 - a(3)) * 0.5);
                q(3) = 0.0;
                q(4) = a(1) * sqrt_tmp;
            end
        end
        
        %for initialization only
        function q = convert_magnetic_field_to_quat(obj, l)
            gamma = l(1)*l(1) + l(2)*l(2);
            sqrt_gamma = sqrt(gamma);
            sqrt_2gamma = sqrt(2.0 * gamma);
            sqrt_2 = sqrt(2.0);

            if l(1) >= 0.0
                sqrt_tmp = sqrt(gamma + l(1)*sqrt_gamma);
                q(1)= sqrt_tmp / sqrt_2gamma;
                q(2) = 0.0;
                q(3) = 0.0;
                q(4) = l(2) / (sqrt_2 * sqrt_tmp);
            else
                sqrt_tmp = sqrt(gamma - l(1)*sqrt_gamma);
                q(1) = l(2) / (sqrt_2 * sqrt_tmp);
                q(2) = 0.0;
                q(3) = 0.0;
                q(4) = sqrt_tmp / sqrt_2gamma;
            end
        end
        
        function ret_obj = init(obj, accel, mag, pz)
            q_accel = convert_gravity_to_quat(obj, -accel);
            q_mag = convert_magnetic_field_to_quat(obj, mag);
            q_init = quaternion_mult(obj, q_accel, q_mag);
            q_init = quaternion_conj(obj, q_init);
            q_init = quat_normalize(obj, q_init);
            obj.x_nominal(7:10) = q_init;
            obj.x_nominal(3) = pz;
            ret_obj = obj;
        end
        
        function ret_obj = predict(obj, ax, ay, az, wx, wy, wz, dt)
            %calculate am - ab
            am_sub_ab = [ax - obj.x_nominal(11);
                         ay - obj.x_nominal(12);
                         az - obj.x_nominal(13)];
            %calculate wm - wb
            wm_sub_wb = [wx - obj.x_nominal(14);
                         wy - obj.x_nominal(15);
                         wz - obj.x_nominal(16)];
            
            %calculate R*(am - ab) + g
            R_am_ab_g = obj.R * am_sub_ab + [0; 0; obj.g_constant];
            
            %nominal state of last time step
            x_last = obj.x_nominal(1:3);
            v_last = obj.x_nominal(4:6);
            q_last = obj.x_nominal(7:10);
            a_b_last = obj.x_nominal(11: 13);
            w_b_last = obj.x_nominal(14: 16);
            
            %first derivative of the quaternion
            q_dot = obj.quaternion_mult(q_last, [0; wm_sub_wb(1); wm_sub_wb(2); wm_sub_wb(3)]);
            
            %nominal state update
            half_dt = 0.5 * dt;
            half_dt_squared = 0.5 * (dt * dt);
            
            %numerical integration
            q_next = q_last + (q_dot .* half_dt);
            q_next = obj.quat_normalize(q_next);
            v_next = v_last + (R_am_ab_g .* dt);
            p_next = x_last + (v_next .* dt) + (R_am_ab_g .* half_dt_squared);
            a_b_next = a_b_last;
            w_b_next = w_b_last;
            
            obj.x_nominal = [p_next(1);
                             p_next(2);
                             p_next(3);
                             v_next(1);
                             v_next(2);
                             v_next(3);
                             q_next(1);
                             q_next(2);
                             q_next(3);
                             q_next(4);
                             a_b_next(1);
                             a_b_next(2);
                             a_b_next(3);
                             w_b_next(1);
                             w_b_next(2);
                             w_b_next(3)];
                         
            %error state update
            F_x = zeros(15, 15);
            F_x(1:3, 1:3) = obj.I_3x3;
            F_x(1:3, 4:6) = obj.I_3x3 .* dt;
            F_x(4:6, 4:6) = obj.I_3x3;
            F_x(4:6, 7:9) = -obj.R * obj.hat_map_3x3(am_sub_ab) * dt;
            F_x(4:6, 10:12) = -obj.R .* dt;
            F_x(7:9, 7:9) = obj.I_3x3 - obj.hat_map_3x3(wm_sub_wb .* dt);
            F_x(7:9, 13:15) = -obj.I_3x3 .* dt;
            F_x(10:12, 10:12) = obj.I_3x3;
            F_x(13:15, 13:15) = obj.I_3x3;
            
            F_i = zeros(15, 12);
            F_i(4:6, 1:3) = obj.I_3x3;
            F_i(7:9, 4:6) = obj.I_3x3;
            F_i(10:12, 7:9) = obj.I_3x3;
            F_i(13:15, 10:12) = obj.I_3x3;
            
            obj.delta_x = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
            %obj.delta_x = F_x * obj.delta_x;
            %disp(obj.delta_x);
            
            obj.P = (F_x * obj.P * F_x.') + (F_i * obj.Q_i * F_i.');
            
            %update rotation matrix for position estimation
            obj.R = obj.quat_to_rotation_matrix(obj.x_nominal(7:10));
            
            ret_obj = obj;
        end
        
        function ret_obj = accel_correct1(obj, ax, ay, az, wx, wy, wz)
            %notice that accelerometer and gyroscope bias are ignored in
            %model1
            
            vx = obj.x_nominal(4);
            vy = obj.x_nominal(5);
            vz = obj.x_nominal(6);
            q0 = obj.x_nominal(7);
            q1 = obj.x_nominal(8);
            q2 = obj.x_nominal(9);
            q3 = obj.x_nominal(10);
            
            %in model1, the measurement vector is gravity instead of
            %accelerometer measurement
            w_cross_v = cross([wx; wy; wz], obj.R.' * [vx; vy; vz]);
            gravity = w_cross_v - [ax; ay; az];
            
            %plot gravity = cross(w, v) - am
            obj.gravity = gravity;
            
            %normalize gravity vector
            y_accel = gravity / norm(gravity);
            
            %error state observation matrix of accelerometer
            H_x_accel = [0 0 0 0 0 0 -2*q2  2*q3 -2*q0 2*q1 0 0 0 0 0 0;
                         0 0 0 0 0 0  2*q1  2*q0  2*q3 2*q2 0 0 0 0 0 0;
                         0 0 0 0 0 0  2*q0 -2*q1 -2*q2 2*q3 0 0 0 0 0 0];
                   
            Q_delte_theta = 0.5 * [-q1 -q2 -q3;
                                    q0 -q3  q2;
                                    q3  q0 -q1;
                                   -q2  q1  q0];
            X_delta_x = zeros(16, 15);
            X_delta_x(1:6, 1:6) = obj.I_6x6;
            X_delta_x(7:10, 7:9) = Q_delte_theta;
            X_delta_x(11:13, 10:12) = obj.I_3x3;
            X_delta_x(14:16, 13:15) = obj.I_3x3;
            
            H_accel = H_x_accel * X_delta_x;

            %prediction of gravity vector using gyroscope
            h_accel = [2*(q1*q3 - q0*q2);
                       2*(q2*q3 + q0*q1);
                       q0*q0 - q1*q1 - q2*q2 + q3*q3];
            
            %calculate kalman gain
            H_accel_t = H_accel.';
            PHt_accel = obj.P * H_accel_t;
            K_accel = PHt_accel / (H_accel * PHt_accel + obj.V_accel);
            %disp(K_accel);
            
            %calculate error state residul
            obj.delta_x = K_accel * (y_accel - h_accel);
            
            %calculate a posteriori process covariance matrix
            obj.P = (obj.I_15x15 - K_accel*H_accel) * obj.P;
            
            %error state injection
            delta_theta_x = obj.delta_x(7);
            delta_theta_y = obj.delta_x(8);
            delta_theta_z = obj.delta_x(9);
            
            if 1
            	q_error = [1;
            	           0.5 * delta_theta_x;
                           0.5 * delta_theta_y;
                           0.5 * delta_theta_z];
            else
                delta_theta_norm = sqrt(delta_theta_x * delta_theta_x + ...
                                        delta_theta_y * delta_theta_y + ...
                                        delta_theta_z * delta_theta_z);
                q_error = [cos(delta_theta_norm / 2);
                           delta_theta_x;
                           delta_theta_y;
                           delta_theta_z];
            end
            
            %error state injection
            %obj.x_nominal(1) = obj.x_nominal(1) + obj.delta_x(1);    %px
            %obj.x_nominal(2) = obj.x_nominal(2) + obj.delta_x(2);    %py
            %obj.x_nominal(3) = obj.x_nominal(3) + obj.delta_x(3);    %pz
            %obj.x_nominal(4) = obj.x_nominal(4) + obj.delta_x(4);    %vx
            %obj.x_nominal(5) = obj.x_nominal(5) + obj.delta_x(5);    %vy
            %obj.x_nominal(6) = obj.x_nominal(6) + obj.delta_x(6);    %vz
            obj.x_nominal(7:10) = obj.quaternion_mult(obj.x_nominal(7:10), q_error); %q
            obj.x_nominal(7:10) = obj.quat_normalize(obj.x_nominal(7:10));           %q
            %obj.x_nominal(11) = obj.x_nominal(11) + obj.delta_x(10); %a_b_x
            %obj.x_nominal(12) = obj.x_nominal(12) + obj.delta_x(11); %a_b_y
            %obj.x_nominal(13) = obj.x_nominal(13) + obj.delta_x(12); %a_b_z
            obj.x_nominal(14) = obj.x_nominal(14) + obj.delta_x(13);  %w_b_x
            obj.x_nominal(15) = obj.x_nominal(15) + obj.delta_x(14);  %w_b_y
            obj.x_nominal(16) = obj.x_nominal(16) + obj.delta_x(15);  %w_b_z
            
            %error state reset
            if 1
                G = obj.I_15x15;
                G(7:9, 7:9) = obj.I_3x3 - (0.5 * obj.hat_map_3x3([delta_theta_x;
                                                                  delta_theta_y;
                                                                  delta_theta_z]));
            else
                G = obj.I_15x15;
            end
            
            obj.P = G * obj.P * G.';
            
            %update rotation matrix for position estimation
            obj.R = obj.quat_to_rotation_matrix(obj.x_nominal(7:10));
            
            ret_obj = obj;
        end
        
        function ret_obj = accel_correct2(obj, ax, ay, az, wx, wy, wz)
            obj.gravity = [0; 0; 0]; %gravity estimation is not valid in acceleromter correction model 2 
            
            vx = obj.x_nominal(4);
            vy = obj.x_nominal(5);
            vz = obj.x_nominal(6);
            q0 = obj.x_nominal(7);
            q1 = obj.x_nominal(8);
            q2 = obj.x_nominal(9);
            q3 = obj.x_nominal(10);
            abx = obj.x_nominal(11);
            aby = obj.x_nominal(12);
            abz = obj.x_nominal(13);
            wbx = obj.x_nominal(14);
            wby = obj.x_nominal(15);
            wbz = obj.x_nominal(16);
            
            %am - ab
            amx_sub_abx = ax - abx;
            amy_sub_aby = ay - aby;
            amz_sub_abz = az - abz;
            
            %wm - wb
            wmx_sub_wbx = wx - wbx;
            wmy_sub_wby = wy - wby;
            wmz_sub_wbz = wz - wbz;
            
            g = obj.g_constant;
            
            %measurement
            y_accel = [amx_sub_abx;
                       amy_sub_aby;
                       amz_sub_abz];
                        
            %error state observation matrix of accelerometer
            dhx_dpx = 0;
            dhy_dpx = 0;
            dhz_dpx = 0;
            dhx_dpy = 0;
            dhy_dpy = 0;
            dhz_dpy = 0;
            dhx_dpz = 0;
            dhy_dpz = 0;
            dhz_dpz = 0;
            dhx_dvx = -wmz_sub_wbz*(2*(q1*q2 - q0*q3)) + wmy_sub_wby*(2*(q1*q3 + q0*q2));
            dhy_dvx = +wmz_sub_wbz*(q0*q0 + q1*q1 - q2*q2 - q3*q3) - wmx_sub_wbx*(2*(q1*q3 + q0*q2));
            dhz_dvx = -wmy_sub_wby*(q0*q0 + q1*q1 - q2*q2 - q3*q3) + wmx_sub_wbx*(2*(q1*q2 - q0*q3));
            dhx_dvy = -wmz_sub_wbz*(q0*q0 - q1*q1 + q2*q2 - q3*q3) + wmy_sub_wby*(2*(q2*q3 - q0*q1));
            dhy_dvy = +wmz_sub_wbz*(2*(q1*q2 + q0*q3)) - wmx_sub_wbx*(2*(q2*q3 - q0*q1));
            dhz_dvy = -wmy_sub_wby*(2*(q1*q2 + q0*q3)) + wmx_sub_wbx*(q0*q0 - q1*q1 + q2*q2 - q3*q3);
            dhx_dvz = -wmz_sub_wbz*(2*(q2*q3+ q0*q1)) + wmy_sub_wby*(q0*q0 - q1*q1 - q2*q2 + q3*q3);
            dhy_dvz = +wmz_sub_wbz*(2*(q1*q3 - q0*q2)) - wmx_sub_wbx*(q0*q0 - q1*q1 - q2*q2 + q3*q3);
            dhz_dvz = -wmy_sub_wby*(2*(q1*q3 - q0*q2)) + wmx_sub_wbx*(2*(q2*q3 + q0*q1));
            dhx_dq0 = -wmz_sub_wbz*(-2*q3*vx + 2*q0*vy + 2*q1*vz) + wmy_sub_wby*(2*q2*vx - 2*q1*vy + 2*q0*vz) + 2*g*q2;
            dhy_dq0 = +wmz_sub_wbz*(2*q0*vx + 2*q3*vy - 2*q2*vz) - wmx_sub_wbx*(2*q2*vx - 2*q1*vy - 2*q0*vz) - 2*g*q1;
            dhz_dq0 = -wmy_sub_wby*(2*q0*vx + 2*q3*vy - 2*q2*vz) + wmx_sub_wbx*(-2*q3*vx + 2*q0*vy + 2*q1*vz) - 2*g*q0;
            dhx_dq1 = -wmz_sub_wbz*(2*q2*vx - 2*q1*vy + 2*q0*vz) + wmy_sub_wby*(2*q3*vx - 2*q0*vy - 2*q1*vz) - 2*g*q3;
            dhy_dq1 = +wmz_sub_wbz*(2*q1*vx + 2*q2*vy + 2*q3*vz) - wmx_sub_wbx*(2*q3*vx - 2*q0*vy - 2*q1*vz) - 2*g*q0;
            dhz_dq1 = -wmy_sub_wby*(2*q1*vx + 2*q2*vy + 2*q3*vz) + wmx_sub_wbx*(2*q2*vx - 2*q1*vy + 2*q0*vz) + 2*g*q1;
            dhx_dq2 = -wmz_sub_wbz*(2*q1*vx + 2*q2*vy + 2*q3*vz) + wmy_sub_wby*(2*q0*vx + 2*q3*vy - 2*q2*vz) + 2*g*q0;
            dhy_dq2 = +wmz_sub_wbz*(-2*q2*vx + 2*q1*vy - 2*q0*vz) - wmx_sub_wbx*(2*q0*vx + 2*q3*vy - 2*q2*vz) - 2*g*q3;
            dhz_dq2 = -wmy_sub_wby*(-2*q2*vx + 2*q1*vy - 2*q0*vz) + wmx_sub_wbx*(2*q1*vx + 2*q2*vy + 2*q3*vz) + 2*g*q2;
            dhx_dq3 = -wmz_sub_wbz*(-2*q0*vx - 2*q3*vy + 2*q2*vz) + wmy_sub_wby*(2*q1*vx + 2*q2*vy + 2*q3*vz) - 2*g*q1;
            dhy_dq3 = +wmz_sub_wbz*(-2*q3*vx + 2*q0*vy + 2*q1*vz) - wmx_sub_wbx*(2*q1*vx + 2*q2*vy + 2*q3*vz) - 2*g*q2;
            dhz_dq3 = -wmy_sub_wby*(-2*q3*vx + 2*q0*vy - 2*q1*vz) + wmx_sub_wbx*(-2*q0*vx - 2*q3*vy + 2*q2*vz) - 2*g*q3;
            dhx_dabx = 0;
            dhy_dabx = 0;
            dhz_dabx = 0;
            dhx_daby = 0;
            dhy_daby = 0;
            dhz_daby = 0;
            dhx_dabz = 0;
            dhy_dabz = 0;
            dhz_dabz = 0;
            dhx_dwbx = 0;
            dhy_dwbx = +(2*(q1*3 + q0*q2)*vx + 2*(q2*q3 - q0*q1)*vy + (q0*q0 - q1*q1 - q2*q2 + q3*q3)*vz);
            dhz_dwbx = -(2*(q1*q2 - q0*q3)*vx + (q0*q0 - q1*q1 + q2*q2 - q3*q3)*vy + 2*(q2*q3 + q0*q1)*vz);
            dhx_dwby = -(2*(q1*q3 + q0*q2)*vx + 2*(q2*q3 - q0*q1)*vy + (q0*q0 - q1*q1 - q2*q2 + q3*q3)*vz);
            dhy_dwby = 0;
            dhz_dwby = +((q0*q0 - q1*q1 - q2*q2 - q3*q3)*vx + 2*(q1*q2 + q0*q3)*vy + 2*(q1*q3 - q0*q2)*vz);
            dhx_dwbz = +(2*(q1*q2 - q0*q3)*vx + (q0*q0 - q1*q1 + q2*q2 - q3*q3)*vy + 2*(q2*q3 + q0*q1)*vz);
            dhy_dwbz = -((q0*q0 + q1*q1 - q2*q2 - q3*q3)*vx + 2*(q1*q2 + q0*q3)*vy + 2*(q1*q3 - q0*q2)*vz);
            dhz_dwbz = 0;
            
            H_x_accel = [dhx_dpx dhx_dpy dhx_dpz dhx_dvx dhx_dvy dhx_dvz dhx_dq0 dhx_dq1 dhx_dq2 dhx_dq3 dhx_dabx dhx_daby dhx_dabz dhx_dwbx dhx_dwby dhx_dwbz;
                         dhy_dpx dhy_dpy dhy_dpz dhy_dvx dhy_dvy dhy_dvz dhy_dq0 dhy_dq1 dhy_dq2 dhy_dq3 dhy_dabx dhy_daby dhy_dabz dhy_dwbx dhy_dwby dhy_dwbz;
                         dhz_dpx dhz_dpy dhz_dpz dhz_dvx dhz_dvy dhz_dvz dhz_dq0 dhz_dq1 dhz_dq2 dhz_dq3 dhz_dabx dhz_daby dhz_dabz dhz_dwbx dhz_dwby dhz_dwbz];
                   
            Q_delte_theta = 0.5 * [-q1 -q2 -q3;
                                    q0 -q3  q2;
                                    q3  q0 -q1;
                                   -q2  q1  q0];
            X_delta_x = zeros(16, 15);
            X_delta_x(1:6, 1:6) = obj.I_6x6;
            X_delta_x(7:10, 7:9) = Q_delte_theta;
            X_delta_x(11:13, 10:12) = obj.I_3x3;
            X_delta_x(14:16, 13:15) = obj.I_3x3;
            
            H_accel = H_x_accel * X_delta_x;

            %prediction of gravity vector using gyroscope
            h_accel = cross([wmx_sub_wbx; wmy_sub_wby; wmz_sub_wbz], obj.R.' * [vx; vy; vz]) - obj.R.' * [0; 0; g];

            %calculate kalman gain
            H_accel_t = H_accel.';
            PHt_accel = obj.P * H_accel_t;
            K_accel = PHt_accel / (H_accel * PHt_accel + obj.V_accel);
            %disp(K_accel);
            
            %calculate error state residul
            obj.delta_x = K_accel * (y_accel - h_accel);
            
            %calculate a posteriori process covariance matrix
            obj.P = (obj.I_15x15 - K_accel*H_accel) * obj.P;
            
            %error state injection
            delta_theta_x = obj.delta_x(7);
            delta_theta_y = obj.delta_x(8);
            delta_theta_z = obj.delta_x(9);
            
            if 1
            	q_error = [1;
            	           0.5 * delta_theta_x;
                           0.5 * delta_theta_y;
                           0.5 * delta_theta_z];
            else
                delta_theta_norm = sqrt(delta_theta_x * delta_theta_x + ...
                                        delta_theta_y * delta_theta_y + ...
                                        delta_theta_z * delta_theta_z);
                q_error = [cos(delta_theta_norm / 2);
                           delta_theta_x;
                           delta_theta_y;
                           delta_theta_z];
            end
            
            %error state injection
            obj.x_nominal(1) = obj.x_nominal(1) + obj.delta_x(1);    %px
            obj.x_nominal(2) = obj.x_nominal(2) + obj.delta_x(2);    %py
            obj.x_nominal(3) = obj.x_nominal(3) + obj.delta_x(3);    %pz
            obj.x_nominal(4) = obj.x_nominal(4) + obj.delta_x(4);    %vx
            obj.x_nominal(5) = obj.x_nominal(5) + obj.delta_x(5);    %vy
            obj.x_nominal(6) = obj.x_nominal(6) + obj.delta_x(6);    %vz
            obj.x_nominal(7:10) = obj.quaternion_mult(obj.x_nominal(7:10), q_error); %q
            obj.x_nominal(7:10) = obj.quat_normalize(obj.x_nominal(7:10));           %q
            obj.x_nominal(11) = obj.x_nominal(11) + obj.delta_x(10); %a_b_x
            obj.x_nominal(12) = obj.x_nominal(12) + obj.delta_x(11); %a_b_y
            obj.x_nominal(13) = obj.x_nominal(13) + obj.delta_x(12); %a_b_z
            obj.x_nominal(14) = obj.x_nominal(14) + obj.delta_x(13); %w_b_x
            obj.x_nominal(15) = obj.x_nominal(15) + obj.delta_x(14); %w_b_y
            obj.x_nominal(16) = obj.x_nominal(16) + obj.delta_x(15); %w_b_z
            
            %error state reset
            if 1
                G = obj.I_15x15;
                G(7:9, 7:9) = obj.I_3x3 - (0.5 * obj.hat_map_3x3([delta_theta_x;
                                                                  delta_theta_y;
                                                                  delta_theta_z]));
            else
                G = obj.I_15x15;
            end
            
            obj.P = G * obj.P * G.';
            
            %update rotation matrix for position estimation
            obj.R = obj.quat_to_rotation_matrix(obj.x_nominal(7:10));
            
            ret_obj = obj;
        end

        function ret_obj = mag_correct1(obj, mx, my, mz)            
            %normalize magnetic field vector
            mag = [mx; my; mz];
            y_mag = mag / norm(mag);
            
            q0 = obj.x_nominal(7);
            q1 = obj.x_nominal(8);
            q2 = obj.x_nominal(9);
            q3 = obj.x_nominal(10);
            
            
            %error state observation matrix of accelerometer
            H_x_mag = [0 0 0 0 0 0 +2*q0 +2*q1 -2*q2 -2*q3 0 0 0 0 0 0;
                       0 0 0 0 0 0 -2*q3 +2*q2 +2*q1 -2*q0 0 0 0 0 0 0;
                       0 0 0 0 0 0 +2*q2 +2*q3 +2*q0 +2*q1 0 0 0 0 0 0];

            Q_delte_theta = 0.5 * [-q1 -q2 -q3;
                                    q0 -q3  q2;
                                    q3  q0 -q1;
                                   -q2  q1  q0];
            X_delta_x = zeros(16, 15);
            X_delta_x(1:6, 1:6) = obj.I_6x6;
            X_delta_x(7:10, 7:9) = Q_delte_theta;
            X_delta_x(11:13, 10:12) = obj.I_3x3;
            X_delta_x(14:16, 13:15) = obj.I_3x3;
            
            H_mag = H_x_mag * X_delta_x;

            %prediction of magnetic field vector using gyroscope
            h_mag = [q0*q0 + q1*q1 - q2*q2 - q3*q3;
                     2*(q1*q2 - q0*q3);
                     2*(q1*q3 + q0*q2)];
            
            %calculate kalman gain
            H_mag_t = H_mag.';
            PHt_mag = obj.P * H_mag_t;
            K_mag = PHt_mag / (H_mag * PHt_mag + obj.V_mag);
            %disp(K_mag);
            
            %calculate error state residul
            obj.delta_x = K_mag * (y_mag - h_mag);
            
            %calculate a posteriori process covariance matrix
            obj.P = (obj.I_15x15 - K_mag*H_mag) * obj.P;
            
            %error state injection
            delta_theta_x = obj.delta_x(7);
            delta_theta_y = obj.delta_x(8);
            delta_theta_z = obj.delta_x(9);
            
            if 1
            	q_error = [1;
            	           0.5 * delta_theta_x;
                           0.5 * delta_theta_y;
                           0.5 * delta_theta_z];
            else
                delta_theta_norm = sqrt(delta_theta_x * delta_theta_x + ...
                                        delta_theta_y * delta_theta_y + ...
                                        delta_theta_z * delta_theta_z);
                q_error = [cos(delta_theta_norm / 2);
                           delta_theta_x;
                           delta_theta_y;
                           delta_theta_z];
            end
            
            %error state injection
            obj.x_nominal(1) = obj.x_nominal(1) + obj.delta_x(1);     %px
            obj.x_nominal(2) = obj.x_nominal(2) + obj.delta_x(2);     %py
            obj.x_nominal(3) = obj.x_nominal(3) + obj.delta_x(3);     %pz
            obj.x_nominal(4) = obj.x_nominal(4) + obj.delta_x(4);     %vx
            obj.x_nominal(5) = obj.x_nominal(5) + obj.delta_x(5);     %vy
            obj.x_nominal(6) = obj.x_nominal(6) + obj.delta_x(6);     %vz
            obj.x_nominal(7:10) = obj.quaternion_mult(obj.x_nominal(7:10), q_error); %q
            obj.x_nominal(7:10) = obj.quat_normalize(obj.x_nominal(7:10));           %q
            obj.x_nominal(11) = obj.x_nominal(11) + obj.delta_x(10);  %a_b_x
            obj.x_nominal(12) = obj.x_nominal(12) + obj.delta_x(11);  %a_b_y
            obj.x_nominal(13) = obj.x_nominal(13) + obj.delta_x(12);  %a_b_z
            obj.x_nominal(14) = obj.x_nominal(14) + obj.delta_x(13);  %w_b_x
            obj.x_nominal(15) = obj.x_nominal(15) + obj.delta_x(14);  %w_b_y
            obj.x_nominal(16) = obj.x_nominal(16) + obj.delta_x(15);  %w_b_z
            
            %error state reset
            if 1
                G = obj.I_15x15;
                G(7:9, 7:9) = obj.I_3x3 - (0.5 * obj.hat_map_3x3([delta_theta_x;
                                                                  delta_theta_y;
                                                                  delta_theta_z]));
            else
                G = obj.I_15x15;
            end
            obj.P = G * obj.P * G.';
            
            %update rotation matrix for position estimation
            obj.R = obj.quat_to_rotation_matrix(obj.x_nominal(7:10));
            
            ret_obj = obj;
        end
        
        function ret_obj = mag_correct2(obj, mx, my, mz)
            %normalize magnetic field vector
            mag = [mx; my; mz];
            y_mag = mag / norm(mag);
            
            q0 = obj.x_nominal(7);
            q1 = obj.x_nominal(8);
            q2 = obj.x_nominal(9);
            q3 = obj.x_nominal(10);
            
            gamma = sqrt(mag(1)*mag(1) + mag(2)*mag(2));
            
            %error state observation matrix of accelerometer
            H_x_mag = [0 0 0 0 0 0 2*(+gamma*q0 - mz*q2) 2*(+gamma*q1 + mz*q3) 2*(-gamma*q2 - mz*q0) 2*(-gamma*q3 + mz*q1) 0 0 0 0 0 0;
                       0 0 0 0 0 0 2*(-gamma*q3 + mz*q1) 2*(+gamma*q2 + mz*q0) 2*(+gamma*q1 + mz*q3) 2*(-gamma*q0 + mz*q2) 0 0 0 0 0 0;
                       0 0 0 0 0 0 2*(+gamma*q2 + mz*q0) 2*(+gamma*q3 - mz*q1) 2*(+gamma*q0 - mz*q2) 2*(+gamma*q1 + mz*q3) 0 0 0 0 0 0];

            Q_delte_theta = 0.5 * [-q1 -q2 -q3;
                                    q0 -q3  q2;
                                    q3  q0 -q1;
                                   -q2  q1  q0];
            X_delta_x = zeros(16, 15);
            X_delta_x(1:6, 1:6) = obj.I_6x6;
            X_delta_x(7:10, 7:9) = Q_delte_theta;
            X_delta_x(11:13, 10:12) = obj.I_3x3;
            X_delta_x(14:16, 13:15) = obj.I_3x3;
            
            H_mag = H_x_mag * X_delta_x;

            %prediction of magnetic field vector using gyroscope
            h_mag = [gamma*(q0*q0 + q1*q1 - q2*q2 - q3*q3) + 2*mz*(q1*q3 - q0*q2);
                     2*(gamma*(q1*q2 - q0*q3) + mz*(q2*q3 + q0*q1));
                     2*gamma*(q1*q3 + q0*q2) + mz*(q0*q0 - q1*q1 - q2*q2 + q3*q3)];
            
            %calculate kalman gain
            H_mag_t = H_mag.';
            PHt_mag = obj.P * H_mag_t;
            K_mag = PHt_mag / (H_mag * PHt_mag + obj.V_mag);
            %disp(K_mag);
            
            %calculate error state residul
            obj.delta_x = K_mag * (y_mag - h_mag);
            
            %calculate a posteriori process covariance matrix
            obj.P = (obj.I_15x15 - K_mag*H_mag) * obj.P;
            
            %error state injection
            delta_theta_x = obj.delta_x(7);
            delta_theta_y = obj.delta_x(8);
            delta_theta_z = obj.delta_x(9);
            
            if 1
            	q_error = [1;
            	           0;
                           0;
                           0.5 * delta_theta_z];
            else
                delta_theta_norm = sqrt(delta_theta_x * delta_theta_x + ...
                                        delta_theta_y * delta_theta_y + ...
                                        delta_theta_z * delta_theta_z);
                q_error = [cos(delta_theta_norm / 2);
                           delta_theta_x;
                           delta_theta_y;
                           delta_theta_z];
            end
            
            %error state injection
            obj.x_nominal(1) = obj.x_nominal(1) + obj.delta_x(1);     %px
            obj.x_nominal(2) = obj.x_nominal(2) + obj.delta_x(2);     %py
            obj.x_nominal(3) = obj.x_nominal(3) + obj.delta_x(3);     %pz
            obj.x_nominal(4) = obj.x_nominal(4) + obj.delta_x(4);     %vx
            obj.x_nominal(5) = obj.x_nominal(5) + obj.delta_x(5);     %vy
            obj.x_nominal(6) = obj.x_nominal(6) + obj.delta_x(6);     %vz
            obj.x_nominal(7:10) = obj.quaternion_mult(obj.x_nominal(7:10), q_error); %q
            obj.x_nominal(7:10) = obj.quat_normalize(obj.x_nominal(7:10));           %q
            obj.x_nominal(11) = obj.x_nominal(11) + obj.delta_x(10);  %a_b_x
            obj.x_nominal(12) = obj.x_nominal(12) + obj.delta_x(11);  %a_b_y
            obj.x_nominal(13) = obj.x_nominal(13) + obj.delta_x(12);  %a_b_z
            obj.x_nominal(14) = obj.x_nominal(14) + obj.delta_x(13);  %w_b_x
            obj.x_nominal(15) = obj.x_nominal(15) + obj.delta_x(14);  %w_b_y
            obj.x_nominal(16) = obj.x_nominal(16) + obj.delta_x(15);  %w_b_z
            
            %error state reset
            if 1
                G = obj.I_15x15;
                G(7:9, 7:9) = obj.I_3x3 - (0.5 * obj.hat_map_3x3([delta_theta_x;
                                                                  delta_theta_y;
                                                                  delta_theta_z]));
            else
                G = obj.I_15x15;
            end
            obj.P = G * obj.P * G.';
            
            %update rotation matrix for position estimation
            obj.R = obj.quat_to_rotation_matrix(obj.x_nominal(7:10));
            
            ret_obj = obj;
        end
        
        function ret_obj = gps_correct(obj, longitude, latitude, vx_ned, vy_ned)
            q0 = obj.x_nominal(7);
            q1 = obj.x_nominal(8);
            q2 = obj.x_nominal(9);
            q3 = obj.x_nominal(10);
            
            %convert longitute/ latitude to ENU position
            pos_ned = obj.covert_geographic_to_ned_frame(longitude, latitude, 0);
            
            %observation vector of gps receiver
            y_gps = [pos_ned(1);
                     pos_ned(2);
                     vx_ned;
                     vy_ned];
                    
            %error state observation matrix of height sensor        
            H_x_gps = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                       0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                       0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
                       0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0];
                      
            Q_delte_theta = 0.5 * [-q1 -q2 -q3;
                                    q0 -q3  q2;
                                    q3  q0 -q1;
                                   -q2  q1  q0];
            X_delta_x = zeros(16, 15);
            X_delta_x(1:6, 1:6) = obj.I_6x6;
            X_delta_x(7:10, 7:9) = Q_delte_theta;
            X_delta_x(11:13, 10:12) = obj.I_3x3;
            X_delta_x(14:16, 13:15) = obj.I_3x3;
            
            H_gps = H_x_gps * X_delta_x;
            
            %prediction of observation vector
            h_gps = H_x_gps * obj.x_nominal;

            %calculate kalman gain
            H_gps_t = H_gps.';
            PHt_gps = obj.P * H_gps_t;
            K_gps = PHt_gps / (H_gps * PHt_gps + obj.V_gps);
            %disp(K_gps);
            
            %calculate error state residul
            obj.delta_x = K_gps * (y_gps - h_gps);
            
            %calculate a posteriori process covariance matrix
            obj.P = (obj.I_15x15 - K_gps*H_gps) * obj.P;
            
            %error state injection
            delta_theta_x = obj.delta_x(7);
            delta_theta_y = obj.delta_x(8);
            delta_theta_z = obj.delta_x(9);
            
            if 1
            	q_error = [1;
            	           0.5 * delta_theta_x;
                           0.5 * delta_theta_y;
                           0.5 * delta_theta_z];
            else
                delta_theta_norm = sqrt(delta_theta_x * delta_theta_x + ...
                                        delta_theta_y * delta_theta_y + ...
                                        delta_theta_z * delta_theta_z);
                q_error = [cos(delta_theta_norm / 2);
                           delta_theta_x;
                           delta_theta_y;
                           delta_theta_z];
            end
            
            %error state injection
            obj.x_nominal(1) = obj.x_nominal(1) + obj.delta_x(1);    %px
            obj.x_nominal(2) = obj.x_nominal(2) + obj.delta_x(2);    %py
            obj.x_nominal(3) = obj.x_nominal(3) + obj.delta_x(3);    %pz
            obj.x_nominal(4) = obj.x_nominal(4) + obj.delta_x(4);    %vx
            obj.x_nominal(5) = obj.x_nominal(5) + obj.delta_x(5);    %vy
            obj.x_nominal(6) = obj.x_nominal(6) + obj.delta_x(6);    %vz
            obj.x_nominal(7:10) = obj.quaternion_mult(obj.x_nominal(7:10), q_error); %q
            obj.x_nominal(7:10) = obj.quat_normalize(obj.x_nominal(7:10));           %q
            obj.x_nominal(11) = obj.x_nominal(11) + obj.delta_x(10); %a_b_x
            obj.x_nominal(12) = obj.x_nominal(12) + obj.delta_x(11); %a_b_y
            obj.x_nominal(13) = obj.x_nominal(13) + obj.delta_x(12); %a_b_z
            obj.x_nominal(14) = obj.x_nominal(14) + obj.delta_x(13); %w_b_x
            obj.x_nominal(15) = obj.x_nominal(15) + obj.delta_x(14); %w_b_y
            obj.x_nominal(16) = obj.x_nominal(16) + obj.delta_x(15); %w_b_z
            
            %error state reset
            G = obj.I_15x15;
            obj.P = G * obj.P * G.';
            
            ret_obj = obj;
        end
        
        function ret_obj = height_correct(obj, pz, vz)
            q0 = obj.x_nominal(7);
            q1 = obj.x_nominal(8);
            q2 = obj.x_nominal(9);
            q3 = obj.x_nominal(10);
            
            %observation vector of height sensor
            y_height = [-pz;
                        -vz];
                    
            %error state observation matrix of height sensor        
            H_x_height = [0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
                          0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0];
                      
            Q_delte_theta = 0.5 * [-q1 -q2 -q3;
                                    q0 -q3  q2;
                                    q3  q0 -q1;
                                   -q2  q1  q0];
            X_delta_x = zeros(16, 15);
            X_delta_x(1:6, 1:6) = obj.I_6x6;
            X_delta_x(7:10, 7:9) = Q_delte_theta;
            X_delta_x(11:13, 10:12) = obj.I_3x3;
            X_delta_x(14:16, 13:15) = obj.I_3x3;
            
            H_height = H_x_height * X_delta_x;
            
            %prediction of observation vector
            h_height = H_x_height * obj.x_nominal;

            %calculate kalman gain
            H_height_t = H_height.';
            PHt_height = obj.P * H_height_t;
            K_height = PHt_height / (H_height * PHt_height + obj.V_height);
            %disp(K_height);
            
            %calculate error state residul
            obj.delta_x = K_height * (y_height - h_height);
            
            %calculate a posteriori process covariance matrix
            obj.P = (obj.I_15x15 - K_height*H_height) * obj.P;
            
            %error state injection
            delta_theta_x = obj.delta_x(7);
            delta_theta_y = obj.delta_x(8);
            delta_theta_z = obj.delta_x(9);
            
            if 1
            	q_error = [1;
            	           0.5 * delta_theta_x;
                           0.5 * delta_theta_y;
                           0.5 * delta_theta_z];
            else
                delta_theta_norm = sqrt(delta_theta_x * delta_theta_x + ...
                                        delta_theta_y * delta_theta_y + ...
                                        delta_theta_z * delta_theta_z);
                q_error = [cos(delta_theta_norm / 2);
                           delta_theta_x;
                           delta_theta_y;
                           delta_theta_z];
            end
            
            %error state injection
            obj.x_nominal(1) = obj.x_nominal(1) + obj.delta_x(1);    %px
            obj.x_nominal(2) = obj.x_nominal(2) + obj.delta_x(2);    %py
            obj.x_nominal(3) = obj.x_nominal(3) + obj.delta_x(3);    %pz
            obj.x_nominal(4) = obj.x_nominal(4) + obj.delta_x(4);    %vx
            obj.x_nominal(5) = obj.x_nominal(5) + obj.delta_x(5);    %vy
            obj.x_nominal(6) = obj.x_nominal(6) + obj.delta_x(6);    %vz
            obj.x_nominal(7:10) = obj.quaternion_mult(obj.x_nominal(7:10), q_error); %q
            obj.x_nominal(7:10) = obj.quat_normalize(obj.x_nominal(7:10));           %q
            obj.x_nominal(11) = obj.x_nominal(11) + obj.delta_x(10); %a_b_x
            obj.x_nominal(12) = obj.x_nominal(12) + obj.delta_x(11); %a_b_y
            obj.x_nominal(13) = obj.x_nominal(13) + obj.delta_x(12); %a_b_z
            obj.x_nominal(14) = obj.x_nominal(14) + obj.delta_x(13); %w_b_x
            obj.x_nominal(15) = obj.x_nominal(15) + obj.delta_x(14); %w_b_y
            obj.x_nominal(16) = obj.x_nominal(16) + obj.delta_x(15); %w_b_z
            
            %error state reset
            G = obj.I_15x15;
            obj.P = G * obj.P * G.';
            
            ret_obj = obj;
        end
    end
end