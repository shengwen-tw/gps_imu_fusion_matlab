classdef eskf
    %assumption: NED coordinate system
    
    properties
        %update time
        dt = 0.001 %100Hz, 0.001s
        
        %nomnial states (16x1)
        p = [0; 0; 0];    %position
        v = [0; 0; 0];    %velocity
        q = [1; 0; 0; 0]; %quaternion
        ab = [0; 0; 0];   %acceleromter bias
        wb = [0; 0; 0];   %gyroscope bias
        g = [0; 0; 9.81]  %gravity vector
        
        %error states (15x1)
        delta_p = [0; 0; 0];        %error position
        delta_v = [0; 0; 0];        %error velocity
        delta_theta = [1; 0; 0; 0]; %error angle vector
        delta_ab = [0; 0; 0];       %error acceleromter bias
        delta_wb = [0; 0; 0];       %error gyroscope bias
        
        %attitude direction cosine matrix
        R = eye(3);  %body-fixed frame to inertial frame
        Rt = eye(3); %inertial frame to body-fixed frame
        
        %noise parameters
        V_i = [0; 0; 0];     %white noise standard deviation of the acceleromter
        Theta_i = [0; 0; 0]; %white noise standard deviation of the gyroscope
        A_i = [0; 0; 0];     %white noise standard deviation of the acceleromater bias
        Omega_i = [0; 0; 0]; %white noise standard deviation of the gyroscope bias
        
        %covariance matrices
        P = eye(16);
        F_i;
        Q_i = eye(16);
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
        
        function nominal_state_update(obj, accelerometer, gyroscope)
            am = accelerometer;
            wm = gyroscope;
            
            %update p
            obj.p = obj.p + (obj.v * dt) + (1/2*(obj.R*(am - obj.ab) + obj.g) * dt * dt);
            
            %update v
            obj.v = v + ((R*(am - obj.ab) + obj.g) * dt);
            
            %update q
            qchange = obj.dt * ...
                [0;
                wm(1) - obj.wb(1);
                wm(2) - obj.wb(2);
                wm(3) - obj.wb(3)];
            q = obj.quaternion_mult(obj.q, q_change);
        end
        
        function error_state_update(obj, accelerometer)
            %update delta_p
            obj.delta_p = obj.delta_p + (obj.delta_v * obj.dt);
            
            %update delta_v
            obj.delta_v = obj.delta_v + (-R*hat_map_3x3(am - obj.ab)*obj.delta_theta - R*obj.delta_ab) * obj.dt;
            
            %update delta_theta
            obj.theta = (obj.Rt * ((obj.wm - obj.wb) * obj.dt) * obj.delta_theta) - ...
                (obj.delta_wb * obj.dt);
            
            %update error state covariance matrix P
        end
        
        function eskf_update(obj, accelerometer)
            obj.nominal_state_update(accelerometer)
            obj.error_state_update(accelerometer)
        end
        
        function nominal_state_injection(obj)
            %inject delta_p into p
            obj.p = obj.p + obj.delta_p;
            
            %inject delta_v into v
            obj.v = obj.v + obj.delta_v;
            
            %inject delta_theta into q
            delta_q = [0; obj.delta_theta(1); obj.delta_theta(2); obj.delta_theta(3)];
            obj.q = obj.quaternion_mult(obj.q, delta_q);
            
            %calculate direction cosine matrix R and Rt from quaternion q
            
            %inject delta_ab_hat into ab
            
            %inject delta_wb_hat into wb
        end
        
        function error_state_reset(obj)
            %reset error state delta_x
            obj.delta_p = [0; 0; 0];
            obj.delta_v = [0; 0; 0];
            obj.delta_theta = [0; 0; 0];
            obj.delta_ab = [0; 0; 0];
            obj.delta_wb = [0; 0; 0];
            
            %reset covariance matrix P
        end
        
        function eskf_correct(obj)
            obj.nominal_state_injection()
            obj.error_state_reset()
        end
    end
end