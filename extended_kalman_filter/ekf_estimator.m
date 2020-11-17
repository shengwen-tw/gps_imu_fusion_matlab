classdef ekf_estimator
    properties
        inclination_angle = 0;
        q_inclination = [1; 0; 0; 0];
        
        %save for last optimal estimation
        x_last = [1; 0; 0; 0]
        
        %prediction a priori state
        x_a_priori = [1; 0; 0; 0]
        %corrected a posterior state
        x_a_posterior = [1; 0; 0; 0]
        
        %process covariance matrix
        P = [2 0 0 0;
             0 2 0 0;
             0 0 2 0;
             0 0 0 2];
        
        %prediction covariance matrix
        Q = [1e-6     0     0     0;
             0     1e-6     0     0;
             0        0  1e-6     0;
             0        0     0  1e-6];
         
       %observation covariance matrix of acceleromter
        R_accel = [1 0 0
                   0 1 0;
                   0 0 1];
        
        %observation covariance matrix of magnetometer
        R_mag = [0.01     0     0;
                 0     0.01     0;
                 0        0  0.01];

        %rotation matrix of current attitude
        R = [1 0 0;
             0 1 0;
             0 0 1];
        
        %euler angles
        roll = 0
        pitch = 0
        yaw = 0
        
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
        
        function ret_obj = predict(obj, wx, wy, wz,dt)
            %update: quaternion ingegration
            half_dt = -0.5 * dt;
            w = [0; wx; wy; wz];
            q_dot = obj.quaternion_mult(w, obj.x_last);
            obj.x_a_priori = [obj.x_last(1) + q_dot(1) * half_dt;
                              obj.x_last(2) + q_dot(2) * half_dt;
                              obj.x_last(3) + q_dot(3) * half_dt;
                              obj.x_last(4) + q_dot(4) * half_dt];
            obj.x_a_priori = obj.quat_normalize(obj.x_a_priori);
            
            F = obj.I_4x4 + (0.5 * dt *...
                [0   -wx  -wy  -wz;
                 wx    0   wz  -wy;
                 wy  -wz    0   wx;
                 wz   wy  -wx    0]);
            
            obj.P = F*obj.P*F.' + obj.Q;
            
            ret_obj = obj;
        end
        
        function ret_obj = correct(obj, gx, gy, gz)
            %normalize gravity vector
            gravity = [gx; gy; gz];
            gravity = gravity / norm(gravity);
            
            q0 = obj.x_a_priori(1);
            q1 = obj.x_a_priori(2);
            q2 = obj.x_a_priori(3);
            q3 = obj.x_a_priori(4);
            
            %correction: acceleromater
            H_accel = [-2*q2   2*q3  -2*q0  2*q1;
                        2*q1   2*q0   2*q3  2*q2;
                        2*q0  -2*q1  -2*q2  2*q3];
            
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
            epsilon_accel(4) = 0;
            %epsilon_accel = obj.quat_normalize(epsilon_accel);
            
            %a posterior estimation
            obj.x_a_posterior = obj.x_a_priori + epsilon_accel;
            obj.x_a_posterior = obj.quat_normalize(obj.x_a_posterior);
            obj.P = (obj.I_4x4 - K_accel*H_accel) * obj.P;
            %disp(obj.P);
            
            obj.x_last = obj.x_a_posterior;
            
            %return the conjugated quaternion since we use the opposite convention compared to the paper
            %paper: quaternion of earth frame to body-fixed frame
            %us: quaternion of body-fixed frame to earth frame
            obj.x_a_posterior = obj.quaternion_conj(obj.x_a_posterior);
            obj.x_a_posterior = obj.quaternion_mult(obj.q_inclination, obj.x_a_posterior);

            %update rotation matrix for position estimation
            obj.R = obj.quat_to_rotation_matrix(obj.x_a_posterior);
            
            %update euler angles for visualization
            euler_angles = obj.quat_to_euler(obj.x_a_posterior);
            obj.roll = euler_angles(1);
            obj.pitch = euler_angles(2);
            obj.yaw = euler_angles(3);
            
            ret_obj = obj;
        end
    end
end