classdef ekf_estimator
    properties
        home_longitude = 0;
        home_latitude = 0;
        home_ecef_x = 0;
        home_ecef_y = 0;
        home_ecef_z = 0;
        
        inclination_angle = 0;
        q_inclination = [1; 0; 0; 0];
        
        %prediction a priori state
        x_a_priori = [0; %px
                      0;  %py
                      0;  %pz
                      0;  %vx
                      0;  %vy
                      0;  %vz
                      1;  %q0
                      0;  %q1
                      0;  %q2
                      0]; %q3
                  
        %corrected a posterior state
        x_a_posterior = [0;  %px
                         0;  %py
                         0;  %pz
                         0;  %vx
                         0;  %vy
                         0;  %vz
                         1;  %q0
                         0;  %q1
                         0;  %q2
                         0]; %q3
        
        %process covariance matrix
        P = [5 0 0 0 0 0 0 0 0 0;  %px
             0 5 0 0 0 0 0 0 0 0;  %py
             0 0 5 0 0 0 0 0 0 0;  %pz
             0 0 0 5 0 0 0 0 0 0;  %vx
             0 0 0 0 5 0 0 0 0 0;  %vy
             0 0 0 0 0 5 0 0 0 0;  %vz
             0 0 0 0 0 0 3 0 0 0;  %q0
             0 0 0 0 0 0 0 3 0 0;  %q1
             0 0 0 0 0 0 0 0 3 0;  %q2
             0 0 0 0 0 0 0 0 0 3]; %q3
        
        %prediction covariance matrix
        Q = [1e-4 0 0 0 0 0 0 0 0 0;  %px
             0 1e-4 0 0 0 0 0 0 0 0;  %py
             0 0 1e-4 0 0 0 0 0 0 0;  %pz
             0 0 0 1e-4 0 0 0 0 0 0;  %vx
             0 0 0 0 1e-4 0 0 0 0 0;  %vy
             0 0 0 0 0 1e-4 0 0 0 0;  %vz
             0 0 0 0 0 0 1e-6 0 0 0;  %q0
             0 0 0 0 0 0 0 1e-6 0 0;  %q1
             0 0 0 0 0 0 0 0 1e-6 0;  %q2
             0 0 0 0 0 0 0 0 0 1e-6]; %q3
         
       %observation covariance matrix of the acceleromter
        R_accel = [0.5 0 0
                   0 0.5 0;
                   0 0 0.5];
        
        %observation covariance matrix of the magnetometer
        R_mag = [10 0 0;
                 0 10 0;
                 0 0 10];

        %observation covariance matrix of the gps sensor
        R_gps = [0.75e-2 0 0 0;  %px
                 0 0.75e-2 0 0;  %py
                 0 0 1.2e-2 0;   %vx
                 0 0 0 1.2e-2];  %vy
        
        %%observation covariance matrix of the height sensor
        R_height = [1e-2 0;  %pz
                    0 1e-2]; %vz
        
        %rotation matrix of current attitude
        R = [1 0 0;
             0 1 0;
             0 0 1];
        
        %euler angles
        roll = 0
        pitch = 0
        yaw = 0
        
        I_3x3 = eye(3);
        I_4x4 = eye(4);
        I_10x10 = eye(10);
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
        
        function ret_obj = set_inclination_angle(obj, angle_degree)
            obj.q_inclination = [cos(deg2rad(angle_degree * 0.5));
                                 0;
                                 0;
                                 sin(deg2rad(angle_degree * 0.5))];
            ret_obj = obj;
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
        
        function quaternion = get_quaternion(obj)
            quaternion = obj.x_a_posterior(7:10);
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

            enu_pos = R * [-dx; -dy; -dz];
        end
        
        function ret_obj = predict(obj, ax, ay, az, wx, wy, wz, dt)
            
            %convert accelerometer's reading from body-fixed frame to
            %inertial frame
            a_inertial = obj.R.' * [ax; ay; az];
            
            %get translational acceleration from accelerometer
            a = [-a_inertial(2);
                 -a_inertial(1);
                 -(a_inertial(3) + 9.8)];
             
            x_last = obj.x_a_posterior(1:3);
            v_last = obj.x_a_posterior(4:6);
            q_last = obj.x_a_posterior(7:10);
            
            %first derivative of the quaternion
            w = [0; wx; wy; wz];
            q_dot = obj.quaternion_mult(w, q_last);

            %state variable update
            half_dt = -0.5 * dt;
            half_dt_squared = 0.5 * (dt * dt);
            q_integration = [q_last(1) + q_dot(1) * half_dt;
                             q_last(2) + q_dot(2) * half_dt;
                             q_last(3) + q_dot(3) * half_dt;
                             q_last(4) + q_dot(4) * half_dt];
            q_integration = obj.quat_normalize(q_integration);
            obj.x_a_priori = [x_last(1) + (v_last(1) * dt) + (a(1) * half_dt_squared); %px
                              x_last(2) + (v_last(2) * dt) + (a(2) * half_dt_squared); %py
                              x_last(3) + (v_last(3) * dt) + (a(3) * half_dt_squared); %pz
                              v_last(1) + (ax * dt); %vx
                              v_last(2) + (ay * dt); %vy
                              v_last(3) + (az * dt); %vz
                              q_integration(1);      %q0
                              q_integration(2);      %q1
                              q_integration(3);      %q2
                              q_integration(4)];     %q3
            
            Omega = obj.I_4x4 + (0.5 * dt *...
                    [0   -wx  -wy  -wz;
                     wx    0   wz  -wy;
                     wy  -wz    0   wx;
                     wz   wy  -wx    0]);
            
            F(1:3, 1:3) = obj.I_3x3;
            F(4:6, 1:3) = obj.I_3x3 .* dt;
            F(4:6, 4:6) = obj.I_3x3;
            F(7:10, 7:10) = Omega;
            
            obj.P = F*obj.P*F.' + obj.Q;
            
            ret_obj = obj;
        end
        
        function ret_obj = accel_correct(obj, gx, gy, gz)
            %normalize gravity vector
            gravity = [gx; gy; gz];
            gravity = gravity / norm(gravity);
            
            q0 = obj.x_a_priori(7);
            q1 = obj.x_a_priori(8);
            q2 = obj.x_a_priori(9);
            q3 = obj.x_a_priori(10);
            
            %correction: acceleromater
            H_accel = [0  0  0  0  0  0  -2*q2   2*q3  -2*q0  2*q1;
                       0  0  0  0  0  0   2*q1   2*q0   2*q3  2*q2;
                       0  0  0  0  0  0   2*q0  -2*q1  -2*q2  2*q3];
            
            %calculate kalman gain
            H_accel_t = H_accel.';
            PHt_accel = obj.P * H_accel_t;
            K_accel = PHt_accel * inv(H_accel * PHt_accel + obj.R_accel);
            %disp(K_accel);
            
            %prediction of gravity vector using gyroscope
            h_accel = [2 * (q1*q3 - q0*q2);
                       2 * (q0*q1 + q2*q3);
                       q0*q0 - q1*q1 - q2*q2 + q3*q3];
                   
            %residual
            epsilon_accel = K_accel * (gravity - h_accel);
            epsilon_accel(10) = 0;  %q3
            epsilon_accel(1:6) = 0; %px, py, pz, vx, vy, vz
            
            %a posterior estimation
            obj.x_a_posterior = obj.x_a_priori + epsilon_accel;
            obj.x_a_posterior(7:10) = obj.quat_normalize(obj.x_a_posterior(7:10));
            obj.P = (obj.I_10x10 - K_accel*H_accel) * obj.P;
            %disp(obj.P);
                        
            %x_a_posterior becomes the x_a_priori of other sensors
            obj.x_a_priori = obj.x_a_posterior;
            
            %update rotation matrix for position estimation
            obj.R = obj.quat_to_rotation_matrix(obj.x_a_posterior(7:10));
            
            %return the conjugated quaternion since we use the opposite convention compared to the paper
            %paper: quaternion of earth frame to body-fixed frame
            %us: quaternion of body-fixed frame to earth frame
            obj.x_a_posterior(7:10) = obj.quaternion_conj(obj.x_a_posterior(7:10));
            
            ret_obj = obj;
        end
        
        function ret_obj = mag_correct(obj, mx, my, mz)
            %normalize gravity vector
            mag = [mx; my; mz];
            mag = mag / norm(mag);
            
            q0 = obj.x_a_priori(7);
            q1 = obj.x_a_priori(8);
            q2 = obj.x_a_priori(9);
            q3 = obj.x_a_priori(10);
            
            %correction: acceleromater
            H_mag = [0  0  0  0  0  0   2*q0  2*q1  -2*q2  -2*q3;
                     0  0  0  0  0  0  -2*q3  2*q2   2*q1  -2*q0;
                     0  0  0  0  0  0   2*q2  2*q3   2*q0   2*q1];
            
            %calculate kalman gain
            H_mag_t = H_mag.';
            PHt_mag = obj.P * H_mag_t;
            K_mag = PHt_mag * inv(H_mag * PHt_mag + obj.R_mag);
            %disp(K_mag);
            
            %prediction of magnetic field vector using gyroscope
            h_mag = [q0*q0 + q1*q1 - q2*q2 - q3*q3;
                     2 * (q1*q2 - q0*q3);
                     2 * (q0*q2 + q1*q3)];
                   
            %residual
            epsilon_mag = K_mag * (mag - h_mag);
            epsilon_mag(8) = 0;   %q2
            epsilon_mag(9) = 0;   %q3
            epsilon_mag(1:6) = 0; %px, py, pz, vx, vy, vz
            
            %a posterior estimation
            obj.x_a_posterior = obj.x_a_priori + epsilon_mag;
            obj.x_a_posterior(7:10) = obj.quat_normalize(obj.x_a_posterior(7:10));
            obj.P = (obj.I_10x10 - K_mag*H_mag) * obj.P;
            %disp(obj.P);
                        
            %x_a_posterior becomes the x_a_priori of other sensors
            obj.x_a_priori = obj.x_a_posterior;
            
            %update rotation matrix for position estimation
            obj.R = obj.quat_to_rotation_matrix(obj.x_a_posterior(7:10));
            
            %return the conjugated quaternion since we use the opposite convention compared to the paper
            %paper: quaternion of earth frame to body-fixed frame
            %us: quaternion of body-fixed frame to earth frame
            obj.x_a_posterior(7:10) = obj.quaternion_conj(obj.x_a_posterior(7:10));
            %obj.x_a_posterior = obj.quaternion_mult(obj.q_inclination, obj.x_a_posterior);
            
            ret_obj = obj;
        end
        
        function ret_obj = gps_correct(obj, longitude, latitude, vx_ned, vy_ned)
            %convert longitute/ latitude to ENU position
            pos_enu = obj.convert_gps_ellipsoid_coordinates_to_enu(longitude, latitude, 0);
            px_enu = pos_enu(1);
            py_enu = pos_enu(2);
            
            %convert NED velocity to ENU frame
            vx_enu = vy_ned;
            vy_enu = vx_ned;
            
            %correction: gps
            H_gps = [1 0 0 0 0 0 0 0 0 0;
                     0 1 0 0 0 0 0 0 0 0;
                     0 0 0 1 0 0 0 0 0 0;
                     0 0 0 0 1 0 0 0 0 0];

            %prediction of observation vector
            h_gps = H_gps * obj.x_a_priori;
                 
            %observation vector of gps sensor
            y_gps = [px_enu;
                     py_enu;
                     vx_enu;
                     vy_enu];
                 
            %calculate kalman gain
            H_gps_t = H_gps.';
            PHt_gps = obj.P * H_gps_t;
            K_gps = PHt_gps * inv(H_gps * PHt_gps + obj.R_gps);
            %disp(K_gps);

            %residual
            residual = K_gps * (y_gps - h_gps);
            residual(3) = 0;     %pz
            residual(6) = 0;     %vz
            residual(7:10) = 0; %q0, q1, q2, q3
            
            %a posterior estimation
            obj.x_a_posterior = obj.x_a_priori + residual;
            obj.P = (obj.I_10x10 - K_gps*H_gps) * obj.P;
            %disp(obj.P);
            
            %x_a_posterior becomes the x_a_priori of other sensors
            obj.x_a_priori = obj.x_a_posterior;
            
            ret_obj = obj;
        end
        
        function ret_obj = height_correct(obj, pz, vz)
            %correction: height sensor
            H_height = [0 0 1 0 0 0 0 0 0 0;
                        0 0 0 0 0 1 0 0 0 0];
            
            %prediction of observation vector
            h_height = H_height * obj.x_a_priori;
                    
            %observation vector of height sensor
            y_height = [pz;
                        vz];
                    
            %calculate kalman gain
            H_height_t = H_height.';
            PHt_height = obj.P * H_height_t;
            K_height = PHt_height * inv(H_height * PHt_height + obj.R_height);
            %disp(K_height);
            
            %residual
            residual = K_height * (y_height - h_height);
            residual(1:2) = 0;  %px, py
            residual(4:5) = 0;  %vx, vy
            residual(7:10) = 0; %q0, q1, q2, q3

            %a posterior estimation
            obj.x_a_posterior = obj.x_a_priori + residual;
            obj.P = (obj.I_10x10 - K_height*H_height) * obj.P;
            %disp(obj.P);

            %x_a_posterior becomes the x_a_priori of other sensors
            obj.x_a_priori = obj.x_a_posterior;
            
            ret_obj = obj;
        end
    end
end