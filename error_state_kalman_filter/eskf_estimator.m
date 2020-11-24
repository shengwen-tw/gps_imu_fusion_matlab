classdef eskf_estimator
    %assumption: NED coordinate system
    
    properties
        home_longitude = 0;
        home_latitude = 0;
        home_ecef_x = 0;
        home_ecef_y = 0;
        home_ecef_z = 0;
        
        inclination_angle = 0;
        q_inclination = [1; 0; 0; 0];
        
        %nomnial state
        x_nominal = [1;  %q0
                     0;  %q2
                     0;  %q2
                     0]; %q3
                      
        %error state
        delta_x = [0;  %theta_x
                   0;  %theta_y
                   0]; %theta_z
                      
        %attitude direction cosine matrix
        R = eye(3);  %inertial frame to body-fixed frame
        Rt = eye(3); %body-fixed frame to inertial frame
        
        %noise parameters
        V_i = [0; 0; 0];     %white noise standard deviation of the acceleromter
        Theta_i = [0.1; 0.1; 0.1]; %white noise standard deviation of the gyroscope
        
        %white noise covariance
        Q_i = [10 0 0;
               0 10 0;
               0 0 10];
        
        %process covariance matrix of error state
        P = [1e-6 0 0;
             0 1e-6 0;
             0 0 1e-6];
        
        %observation covariance matrix of accelerometer
        V_accel = [10 0 0;
                   0 10 0;
                   0 0 10];
         
        I_3x3 = eye(3);
        I_4x4 = eye(4);
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
            quaternion = obj.x_nominal(1:4);
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

            enu_pos = R * [-dx; -dy; -dz];
        end
        
        function ret_obj = predict(obj, wx, wy, wz, dt)
            q_last = obj.x_nominal(1:4);
            
            %first derivative of the quaternion
            w = [0; wx; wy; wz];
            q_dot = obj.quaternion_mult(w, q_last);
            
            %nominal state update
            half_dt = -0.5 * dt;
            %half_dt_squared = 0.5 * (dt * dt);
            q_integration = [q_last(1) + q_dot(1) * half_dt;
                             q_last(2) + q_dot(2) * half_dt;
                             q_last(3) + q_dot(3) * half_dt;
                             q_last(4) + q_dot(4) * half_dt];
            q_integration = obj.quat_normalize(q_integration);
            obj.x_nominal = [q_integration(1);  %q0
                             q_integration(2);  %q1
                             q_integration(3);  %q2
                             q_integration(4)]; %q3
                         
            %error state update
            F_x = obj.I_3x3 - obj.hat_map_3x3([wx; wy; wz] .* dt);
            F_i = eye(3);
            obj.delta_x = F_x * obj.delta_x;
            obj.P = (F_x * obj.P * F_x.') + (F_i * obj.Q_i * F_i.');
            
            %update rotation matrix for position estimation
            obj.R = obj.quat_to_rotation_matrix(obj.x_nominal(1:4));
            
            ret_obj = obj;
        end
        
        function ret_obj = accel_correct(obj, gx, gy, gz)
            %normalize gravity vector
            gravity = [gx; gy; gz];
            y = gravity / norm(gravity);
            
            q0 = obj.x_nominal(1);
            q1 = obj.x_nominal(2);
            q2 = obj.x_nominal(3);
            q3 = obj.x_nominal(4);
            
            %error state observation matrix of accelerometer
            H_x_accel = [-2*q2   2*q3  -2*q0  2*q1;
                          2*q1   2*q0   2*q3  2*q2;
                          2*q0  -2*q1  -2*q2  2*q3];
                   
            Q_delte_theta = 0.5 * [-q1 -q2 -q3;
                                    q0 -q3  q2;
                                    q3  q0 -q1;
                                   -q2  q1  q0];                
            X_delta_x = Q_delte_theta;
            H_accel = H_x_accel * X_delta_x;

            %prediction of gravity vector using gyroscope
            h_accel = [2 * (q1*q3 - q0*q2);
                       2 * (q0*q1 + q2*q3);
                       q0*q0 - q1*q1 - q2*q2 + q3*q3];
            
            %calculate kalman gain
            H_accel_t = H_accel.';
            PHt_accel = obj.P * H_accel_t;
            K_accel = PHt_accel * inv(H_accel * PHt_accel + obj.V_accel);
            %disp(K_accel);
            
            %calculate error state residul
            obj.delta_x = K_accel * (y - h_accel);
            
            %calculate a posteriori process covariance matrix
            obj.P = (obj.I_3x3 - K_accel*H_accel) * obj.P;
            
            %error state injection
            delta_theta_x = obj.delta_x(1);
            delta_theta_y = obj.delta_x(2);
            delta_theta_z = obj.delta_x(3);
            %delta_theta_norm = sqrt(delta_theta_x * delta_theta_x + ...
            %                        delta_theta_y * delta_theta_y + ...
            %                        delta_theta_z * delta_theta_z);
            %q_error = [cos(delta_theta_norm / 2);
            %           delta_theta_x;
            %           delta_theta_y;
            %           delta_theta_z];
            q_error = [1;
                       0.5 * delta_theta_x;
                       0.5 * delta_theta_y;
                       0.5 * delta_theta_z];
            obj.x_nominal(1:4) = obj.quaternion_mult(obj.x_nominal(1:4), q_error);
            
            %error state reset
            obj.delta_x = [0; 0; 0];
            %G = obj.I_3x3 - (0.5 * hat_map_3x3([delta_theta_x;
            %                                    delta_theta_y;
            %                                    delta_theta_z]));
            G = obj.I_3x3;
            obj.P = G * obj.P * G.';
            
            %update rotation matrix for position estimation
            obj.R = obj.quat_to_rotation_matrix(obj.x_nominal(1:4));
            
            ret_obj = obj;
        end
    end
end