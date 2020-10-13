classdef attitude
    properties
        q_last = [1; 0; 0; 0]
        q_estimate = [1; 0; 0; 0]
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
            norm = sqrt(q(1)*q(1) + q(2)*q(2) + q(3)*q(3) + q(4)*q(4));
            q_out = [q(1);
                     q(2);
                     q(3);
                     q(4)];
        end
        
        function q_out = lerp(obj, q1, q2, alpha)
            q_out = [(1.0 - alpha) * q1(1) + (alpha * q2(1));
                     (1.0 - alpha) * q1(2) + (alpha * q2(2));
                     (1.0 - alpha) * q1(3) + (alpha * q2(3));
                     (1.0 - alpha) * q1(4) + (alpha * q2(4))];
        end
        
        function q = convert_gravity_to_quat(obj, gx, gy, gz)
            if gz >= 0
                q = [sqrt((gx + 1) / 2);
                     -gy / sqrt(2 * (gz + 1));
                     gx / sqrt(2 * (gz + 1));
                     0];
            else
                q = [-gy / sqrt(2 * (1 - gz));
                     sqrt((1 - gz) / 2);
                     0;
                     gx / sqrt(2 * (1 - gz))];
            end
        end
        
        function q = convert_magnetic_field_to_quat(obj, mx, my, mz)
            gamma = (mx * mx) + (my * my);
            sqrt_gamma = sqrt(gamma);
            sqrt_2_gamma = sqrt(2 * gamma);
            sqrt_2 = sqrt(2);
            
            if mx >= 0
                q = [sqrt(gamma + (mx * sqrt_gamma)) / sqrt_2_gamma;
                     0;
                     0;
                     my / (sqrt_2 * sqrt(gamma + (mx * sqrt_gamma)))];
            else
                q = [my / (sqrt_2 * sqrt(gamma - (mx * sqrt_gamma)));
                     0;
                     0;
                     sqrt(gamma - mx * (sqrt_gamma)) / sqrt_2_gamma];
            end
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

        function ret_obj = complementary_filter(obj, gx, gy, gz, wx, wy, wz, mx, my, mz, dt)
            %constant
            q_identy = [1; 0; 0; 0];
            
            %quaternion ingegration
            half_dt = 0.5 * dt;
            w = [0; wx; wy; wz];
            q_dot = obj.quaternion_mult(w, obj.q_last);
            q_gyro = [obj.q_last(1) + q_dot(1) * half_dt;
                      obj.q_last(2) + q_dot(2) * half_dt;
                      obj.q_last(3) + q_dot(3) * half_dt;
                      obj.q_last(4) + q_dot(4) * half_dt];
            
            %calculate predicted gravity vector
            conj_q_gyro = obj.quaternion_conj(q_gyro);
            R_gyro = obj.prepare_body_to_earth_rotation_matrix(conj_q_gyro);
            g_predict = R_gyro * [gx; gy; gz]; 
            g_predict = g_predict / norm(g_predict);
            
            %calculate predicted magnetic field vector
            l_predict = R_gyro * [mx; my; mz]; 
            l_predict = l_predict / norm(l_predict);
            
            %calculate delta change of quaternion for fusing gyroscope and aceelerometer
            weight_accel = 0.005;
            delta_q_acc = obj.convert_gravity_to_quat(g_predict(1), g_predict(2), g_predict(3));
            delta_q_acc = obj.quat_normalize(delta_q_acc);
            bar_delta_q_acc = obj.lerp(q_identy, delta_q_acc, weight_accel);
            bar_delta_q_acc = obj.quat_normalize(bar_delta_q_acc);
            
            %calculate delta change of quaternion for fusing gyroscope and magnetometer
            weight_mag = 0.005;
            delta_q_mag = obj.convert_magnetic_field_to_quat(l_predict(1), l_predict(2), l_predict(3));
            delta_q_mag = obj.quat_normalize(delta_q_mag);
            bar_delta_q_mag = obj.lerp(q_identy, delta_q_mag, weight_mag);
            bar_delta_q_mag = obj.quat_normalize(bar_delta_q_mag);
            
            %fuse gyroscope, accelerometer and magnetometer using delte change
            q_delta_acc_mag = obj.quaternion_mult(bar_delta_q_acc, bar_delta_q_mag);
            obj.q_estimate = obj.quaternion_mult(q_gyro, q_delta_acc_mag);
            obj.q_last = obj.q_estimate;
            
            %return the conjugated quaternion since we use the opposite convention compared to the paper
            %paper: quaternion of earth frame to body-fixed frame
            %us: quaternion of body-fixed frame to earth frame
            obj.q_estimate = obj.quaternion_conj(obj.q_estimate);
            
            ret_obj = obj;
        end
    end
end