classdef eskf_estimator
    %assumption: NED coordinate system
    
    properties
        home_longitude = 0;
        home_latitude = 0;
        home_ecef_x = 0;
        home_ecef_y = 0;
        home_ecef_z = 0;
        
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
                     0]; %q3
                      
        %error state
        delta_x = [0;  %px
                   0;  %py
                   0;  %pz
                   0;  %vx
                   0;  %vy
                   0;  %vz
                   0;  %theta_x
                   0;  %theta_y
                   0]; %theta_z
                      
        %attitude direction cosine matrix
        R = eye(3);  %inertial frame to body-fixed frame
        Rt = eye(3); %body-fixed frame to inertial frame
        
        %noise parameters
        V_i = [];     %white noise standard deviation of the acceleromter
        Theta_i = []; %white noise standard deviation of the gyroscope
        
        %covariance matrix of process white noise
        Q_i = [1e-6 0 0 0 0 0;  %noise of ax
               0 1e-6 0 0 0 0;  %noise of ay
               0 0 1e-6 0 0 0;  %noise of az
               0 0 0 1e-5 0 0;  %noise of wx
               0 0 0 0 1e-5 0;  %noise of wy
               0 0 0 0 0 1e-5]; %noise of wz
        
        %process covariance matrix of error state
        P = [5 0 0 0 0 0 0 0 0;  %px
             0 5 0 0 0 0 0 0 0;  %py
             0 0 5 0 0 0 0 0 0;  %pz
             0 0 0 5 0 0 0 0 0;  %vx
             0 0 0 0 5 0 0 0 0;  %vy
             0 0 0 0 0 5 0 0 0;  %vz
             0 0 0 0 0 0 5 0 0;  %delta_x
             0 0 0 0 0 0 0 5 0;  %delta_y
             0 0 0 0 0 0 0 0 5]; %delta_z
        
        %observation covariance matrix of accelerometer
        V_accel = [7e-2 0 0;  %ax
                   0 7e-2 0;  %ay
                   0 0 7e-2]; %az
               
        %observation covariance matrix of accelerometer
        V_mag = [5 0 0;  %mx
                 0 5 0;  %my
                 0 0 5]; %mz
             
        %observation covariance matrix of the gps sensor
        V_gps = [5e-5 0 0 0;  %px
                 0 5e-5 0 0;  %py
                 0 0 1e-4 0;   %vx
                 0 0 0 1e-4];  %vy
             
        %%observation covariance matrix of the height sensor
        V_height = [1e-5 0;  %pz
                    0 1e-2]; %vz   
                
        I_3x3 = eye(3);
        I_4x4 = eye(4);
        I_6x6 = eye(6);
        I_9x9 = eye(9);
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

        function R = prepare_body_to_earth_rotation_matrix(obj, q)
            q0q0 = q(1) * q(1);
            q1q1 = q(2) * q(2);
            q2q2 = q(3) * q(3);
            q3q3 = q(4) * q(4);
            q0q3 = q(1) * q(4);
            q1q2 = q(2) * q(3);
            q2q3 = q(3) * q(4);
            q1q3 = q(2) * q(4);
            q0q2 = q(1) * q(3);
            q0q1 = q(1) * q(2);

            r11 = q0q0 + q1q1 - q2q2 - q3q3;
            r12 = 2*(q1q2 - q0q3);
            r13 = 2*(q1q3 + q0q2);

            r21 = 2*(q1q2 + q0q3);
            r22 = q0q0 - q1q1 + q2q2 - q3q3;
            r23 = 2*(q2q3 - q0q1);

            r31 = 2*(q1q3 - q0q2);
            r32 = 2*(q2q3 + q0q1);
            r33 = q0q0 - q1q1 - q2q2 + q3q3;
            
            R = [r11 r12 r13;
                 r21 r22 r23;
                 r31 r32 r33];
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
            %return the conjugated quaternion since we use the opposite convention compared to the paper
            %paper: quaternion of earth frame to body-fixed frame
            %us: quaternion of body-fixed frame to earth frame
            quaternion = obj.x_nominal(7:10);
            quaternion = obj.quaternion_conj(quaternion);
            %quaternion = obj.quaternion_mult(obj.q_inclination, quaternion);
        end
        function vec_enu = convert_3x1_vector_ned_to_enu(obj, vec_ned)
            vec_enu = [+vec_ned(2);
                       +vec_ned(1);
                       -vec_ned(3)];
        end
        
        function ret_obj = set_home_longitude_latitude(obj, longitude, latitude, height_msl)
            EARTH_RADIUS = 6371000;
            
            sin_lambda = sin(deg2rad(longitude));
            cos_lambda = cos(deg2rad(longitude));
            sin_phi = sin(deg2rad(latitude));
            cos_phi = cos(deg2rad(latitude));

            obj.home_longitude = longitude;
            obj.home_latitude = latitude;
            obj.home_ecef_x = (height_msl + EARTH_RADIUS) * cos_phi * cos_lambda;
            obj.home_ecef_y = (height_msl + EARTH_RADIUS) * cos_phi * sin_lambda;
            obj.home_ecef_z = (height_msl + EARTH_RADIUS) * sin_phi;
            
            ret_obj = obj;
        end
        
        function enu_pos = convert_gps_ellipsoid_coordinates_to_enu(obj, longitude, latitude, height_msl)
            EARTH_RADIUS = 6371000;
            
            sin_lambda = sin(deg2rad(longitude));
            cos_lambda = cos(deg2rad(longitude));
            sin_phi = sin(deg2rad(latitude));
            cos_phi = cos(deg2rad(latitude));

            %convert geodatic coordinates to earth center earth fixed frame (ecef)
            ecef_now_x = (height_msl + EARTH_RADIUS) * cos_phi * cos_lambda;
            ecef_now_y = (height_msl + EARTH_RADIUS) * cos_phi * sin_lambda;
            ecef_now_z = (height_msl + EARTH_RADIUS) * sin_phi;
            
            %convert position from earth center earth fixed frame to east north up frame
            r11 = -sin_lambda;
            r12 = cos_lambda;
            r13 = 0;
            r21 = -cos_lambda * sin_phi;
            r22 = -sin_lambda * sin_phi;
            r23 = cos_phi;
            r31 = cos_lambda * cos_phi;
            r32 = sin_lambda * cos_phi;
            r33 = sin_phi;
            R = [r11, r12, r13;
                 r21, r22, r23
                 r31, r32, r33];
            
            dx = ecef_now_x - obj.home_ecef_x;
            dy = ecef_now_y - obj.home_ecef_y;
            dz = ecef_now_z - obj.home_ecef_z;

            enu_pos = R * [dx; dy; dz];
        end
        
        function ret_obj = predict(obj, ax, ay, az, wx, wy, wz, dt)
            %convert accelerometer's reading from body-fixed frame to
            %inertial frame
            a_inertial = obj.R.' * [ax; ay; az];
            
            %get translational acceleration from accelerometer
            a_ned = [a_inertial(1);
                     a_inertial(2);
                     a_inertial(3) + 9.8];
            
            a = [a_ned(2); a_ned(1); -a_ned(3)]; %NED to ENU
             
            x_last = obj.x_nominal(1:3);
            v_last = obj.x_nominal(4:6);
            q_last = obj.x_nominal(7:10);
            
            %first derivative of the quaternion
            w = [0; wx; wy; wz];
            q_dot = obj.quaternion_mult(w, q_last);
            
            %nominal state update
            half_dt = -0.5 * dt;
            half_dt_squared = 0.5 * (dt * dt);
            q_integration = [q_last(1) + q_dot(1) * half_dt;
                             q_last(2) + q_dot(2) * half_dt;
                             q_last(3) + q_dot(3) * half_dt;
                             q_last(4) + q_dot(4) * half_dt];
            q_integration = obj.quat_normalize(q_integration);
            obj.x_nominal = [x_last(1) + (v_last(1) * dt) + (a(1) * half_dt_squared); %px
                             x_last(2) + (v_last(2) * dt) + (a(2) * half_dt_squared); %py
                             x_last(3) + (v_last(3) * dt) + (a(3) * half_dt_squared); %pz
                             v_last(1) + (a(1) * dt); %vx
                             v_last(2) + (a(2) * dt); %vy
                             v_last(3) + (a(3) * dt); %vz
                             q_integration(1);      %q0
                             q_integration(2);      %q1
                             q_integration(3);      %q2
                             q_integration(4)];     %q3
                         
            %error state update
            F_x = zeros(9, 9);
            F_x(1:3, 1:3) = obj.I_3x3;
            F_x(1:3, 4:6) = obj.I_3x3 .* dt;
            F_x(4:6, 4:6) = obj.I_3x3;
            F_x(4:6, 7:9) = -obj.R * obj.hat_map_3x3(a) * dt;
            F_x(7:9, 7:9) = obj.I_3x3 - obj.hat_map_3x3([wx; wy; wz] .* dt);
            
            F_i = zeros(9, 6);
            F_i(4:6, 1:3) = obj.I_3x3;
            F_i(7:9, 4:6) = obj.I_3x3;
            
            obj.P = (F_x * obj.P * F_x.') + (F_i * obj.Q_i * F_i.');
            
            %update rotation matrix for position estimation
            obj.R = obj.quat_to_rotation_matrix(obj.x_nominal(7:10));
            
            ret_obj = obj;
        end
        
        function ret_obj = accel_correct(obj, gx, gy, gz)
            %normalize gravity vector
            gravity = [gx; gy; gz];
            y_accel = gravity / norm(gravity);
            
            q0 = obj.x_nominal(7);
            q1 = obj.x_nominal(8);
            q2 = obj.x_nominal(9);
            q3 = obj.x_nominal(10);
            
            %error state observation matrix of accelerometer
            H_x_accel = [0 0 0 0 0 0  2*q2   2*q3   2*q0  2*q1;
                         0 0 0 0 0 0 -2*q1  -2*q0   2*q3  2*q2;
                         0 0 0 0 0 0  2*q0  -2*q1  -2*q2  2*q3];
                   
            Q_delte_theta = 0.5 * [-q1 -q2 -q3;
                                    q0 -q3  q2;
                                    q3  q0 -q1;
                                   -q2  q1  q0];
            X_delta_x = zeros(10, 9);
            X_delta_x(1:6, 1:6) = obj.I_6x6;
            X_delta_x(7:10, 7:9) = Q_delte_theta;
            
            H_accel = H_x_accel * X_delta_x;

            %prediction of gravity vector using gyroscope
            h_accel = [2 * (q0*q2 + q1*q3);
                       2 * (q2*q3 - q0*q1);
                       q0*q0 - q1*q1 - q2*q2 + q3*q3];
            
            %calculate kalman gain
            H_accel_t = H_accel.';
            PHt_accel = obj.P * H_accel_t;
            K_accel = PHt_accel / (H_accel * PHt_accel + obj.V_accel);
            %disp(K_accel);
            
            %calculate error state residul
            obj.delta_x = K_accel * (y_accel - h_accel);
            
            %calculate a posteriori process covariance matrix
            obj.P = (obj.I_9x9 - K_accel*H_accel) * obj.P;
            
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
            
            obj.x_nominal(7:10) = obj.quaternion_mult(obj.x_nominal(7:10), q_error);
            obj.x_nominal(7:10) = obj.quat_normalize(obj.x_nominal(7:10));
            
            %error state reset
            if 1
                G = obj.I_9x9;
                G(7:9, 7:9) = obj.I_3x3 - (0.5 * obj.hat_map_3x3([delta_theta_x;
                                                                  delta_theta_y;
                                                                  delta_theta_z]));
            else
                G = obj.I_9x9;
            end
            
            obj.P = G * obj.P * G.';
            
            %update rotation matrix for position estimation
            obj.R = obj.quat_to_rotation_matrix(obj.x_nominal(7:10));
            
            ret_obj = obj;
        end
        
        function ret_obj = mag_correct(obj, mx, my, mz)
            %normalize magnetic field vector
            mag = [mx; my; mz];
            y_mag = mag / norm(mag);
            
            q0 = obj.x_nominal(7);
            q1 = obj.x_nominal(8);
            q2 = obj.x_nominal(9);
            q3 = obj.x_nominal(10);
            
            gamma = sqrt(mag(1)*mag(1) + mag(2)*mag(2));
            
            %error state observation matrix of accelerometer
            H_x_mag = [0 0 0 0 0 0 2*(+gamma*q0 + mag(3)*q2) 2*(+gamma*q1 + mag(3)*q3) 2*(-gamma*q2 + mag(3)*q0) ...
                                   2*(-gamma*q3 + mag(3)*q1);
                       0 0 0 0 0 0 2*(+gamma*q3 - mag(3)*q1) 2*(+gamma*q2 - mag(3)*q0) 2*(+gamma*q1 + mag(3)*q3) ...
                                   2*(+gamma*q0 + mag(3)*q2);
                       0 0 0 0 0 0 2*(-gamma*q2 + mag(3)*q0) 2*(+gamma*q3 - mag(3)*q1) 2*(-gamma*q0 - mag(3)*q2) ...
                                   2*(+gamma*q1 + mag(3)*q3)];

            Q_delte_theta = 0.5 * [-q1 -q2 -q3;
                                    q0 -q3  q2;
                                    q3  q0 -q1;
                                   -q2  q1  q0];
            X_delta_x = zeros(10, 9);
            X_delta_x(1:6, 1:6) = obj.I_6x6;
            X_delta_x(7:10, 7:9) = Q_delte_theta;

            H_mag = H_x_mag * X_delta_x;

            %prediction of magnetic field vector using gyroscope
            h_mag = [gamma*(q0*q0 + q1*q1 - q2*q2 - q3*q3) + 2*mag(3)*(q0*q2 + q1*q3);
                     2*gamma*(q1*q2 + q0*q3) + 2*mag(3)*(q2*q3 - q0*q1);
                     2*gamma*(q1*q3 - q0*q2) + mag(3)*(q0*q0 - q1*q1 - q2*q2 + q3*q3)];
            
            %calculate kalman gain
            H_mag_t = H_mag.';
            PHt_mag = obj.P * H_mag_t;
            K_mag = PHt_mag / (H_mag * PHt_mag + obj.V_mag);
            %disp(K_mag);
            
            %calculate error state residul
            obj.delta_x = K_mag * (y_mag - h_mag);
            
            %calculate a posteriori process covariance matrix
            obj.P = (obj.I_9x9 - K_mag*H_mag) * obj.P;
            
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
            
            obj.x_nominal(7:10) = obj.quaternion_mult(obj.x_nominal(7:10), q_error);
            obj.x_nominal(7:10) = obj.quat_normalize(obj.x_nominal(7:10));
            
            %error state reset
            if 1
                G = obj.I_9x9;
                G(7:9, 7:9) = obj.I_3x3 - (0.5 * obj.hat_map_3x3([delta_theta_x;
                                                                  delta_theta_y;
                                                                  delta_theta_z]));
            else
                G = obj.I_9x9;
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
            pos_enu = obj.convert_gps_ellipsoid_coordinates_to_enu(longitude, latitude, 0);
            px_enu = pos_enu(1);
            py_enu = pos_enu(2);
            
            %convert NED velocity to ENU frame
            vx_enu = vy_ned;
            vy_enu = vx_ned;
            
            %observation vector of gps receiver
            y_gps = [px_enu;
                     py_enu;
                     vx_enu;
                     vy_enu];
                    
            %error state observation matrix of height sensor        
            H_x_gps = [1 0 0 0 0 0 0 0 0 0;
                       0 1 0 0 0 0 0 0 0 0;
                       0 0 0 1 0 0 0 0 0 0;
                       0 0 0 0 1 0 0 0 0 0];
                      
            Q_delte_theta = 0.5 * [-q1 -q2 -q3;
                                    q0 -q3  q2;
                                    q3  q0 -q1;
                                   -q2  q1  q0];
            X_delta_x = zeros(10, 9);
            X_delta_x(1:6, 1:6) = obj.I_6x6;
            X_delta_x(7:10, 7:9) = Q_delte_theta;

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
            obj.P = (obj.I_9x9 - K_gps*H_gps) * obj.P;
            
            %error state injection
            obj.x_nominal(1) = obj.x_nominal(1) + obj.delta_x(1);
            obj.x_nominal(2) = obj.x_nominal(2) + obj.delta_x(2);
            obj.x_nominal(4) = obj.x_nominal(4) + obj.delta_x(4);
            obj.x_nominal(5) = obj.x_nominal(5) + obj.delta_x(5);
            
            %error state reset
            G = obj.I_9x9;
            obj.P = G * obj.P * G.';
            
            ret_obj = obj;
        end
        
        function ret_obj = height_correct(obj, pz, vz)
            q0 = obj.x_nominal(7);
            q1 = obj.x_nominal(8);
            q2 = obj.x_nominal(9);
            q3 = obj.x_nominal(10);
            
            %observation vector of height sensor
            y_height = [pz;
                        vz];
                    
            %error state observation matrix of height sensor        
            H_x_height = [0 0 1 0 0 0 0 0 0 0;
                          0 0 0 0 0 1 0 0 0 0];
                      
            Q_delte_theta = 0.5 * [-q1 -q2 -q3;
                                    q0 -q3  q2;
                                    q3  q0 -q1;
                                   -q2  q1  q0];
            X_delta_x = zeros(10, 9);
            X_delta_x(1:6, 1:6) = obj.I_6x6;
            X_delta_x(7:10, 7:9) = Q_delte_theta;

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
            obj.P = (obj.I_9x9 - K_height*H_height) * obj.P;
            
            %error state injection
            obj.x_nominal(3) = obj.x_nominal(3) + obj.delta_x(3);
            obj.x_nominal(6) = obj.x_nominal(6) + obj.delta_x(6);
            
            %error state reset
            G = obj.I_9x9;
            obj.P = G * obj.P * G.';
            
            ret_obj = obj;
        end
    end
end