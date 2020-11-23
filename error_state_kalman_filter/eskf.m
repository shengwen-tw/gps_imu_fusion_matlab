classdef eskf
    %assumption: NED coordinate system
    
    properties 
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
        Theta_i = [0; 0; 0]; %white noise standard deviation of the gyroscope
        
        %white noise covariance
        Q_i = [Theta_i(1) 0 0;
               0 Theta_i(2) 0;
               0 0 Theta_i(3)];
        
        %process covariance matrix of error state
        P = [1e-3 0 0;
             0 1e-3 0;
             0 0 1e-3];
        
        %observation covariance matrix of accelerometer
        V_accel = [0.5 0 0;
                   0 0.5 0;
                   0 0 0.5];
         
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
        
        function skew_matrix = hat_map_3x3(skew_vector)
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
        
        function ret_obj = predict(obj, wx, wy, wz)
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
            F_x = obj.Rt * ([wx; wy; wz] .* dt);
            F_i = eye(3);
            obj.delta_x = F_x * obj.delta_x;
            obj.P = (F_x * obj.P * F_x.') + (F_i * obj.Q * F_i.');
            
            ret_obj = obj;
        end
        
        function ret_obj = accelerometer_correct(obj, gx, gy, gz)
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
            %delta_theta_x = obj.delta_x(1);
            %delta_theta_y = obj.delta_x(2);
            %delta_theta_z = obj.delta_x(3);
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
            G = obj.I_4x4 - (0.5 * dt * ...
                    [0   -wx  -wy  -wz;
                     wx    0   wz  -wy;
                     wy  -wz    0   wx;
                     wz   wy  -wx    0]);
            obj.P = G * obj.P * G.';
            
            ret_obj = obj;
        end
    end
end