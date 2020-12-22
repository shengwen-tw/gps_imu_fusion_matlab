classdef attitude_estimator
    properties        
        q_last = [1; 0; 0; 0]
        q_estimate = [1; 0; 0; 0]
        
        R = [1 0 0;
             0 1 0;
             0 0 1];
        
        roll = 0
        pitch = 0
        yaw = 0
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
        
        function q_out = lerp(obj, q1, q2, alpha)
            q_out = [(1.0 - alpha) * q1(1) + (alpha * q2(1));
                     (1.0 - alpha) * q1(2) + (alpha * q2(2));
                     (1.0 - alpha) * q1(3) + (alpha * q2(3));
                     (1.0 - alpha) * q1(4) + (alpha * q2(4))];
        end
        
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
        
        function ret_obj = complementary_filter(obj, gx, gy, gz, wx, wy, wz, mx, my, mz, dt)
            %constant
            q_identy = [1; 0; 0; 0];
            
            %normalize accelerometer and magnetometer's reading
            accelerometer = [gx; gy; gz];
            accelerometer = accelerometer / norm(accelerometer);
            magnetometer = [mx; my; mz];
            magnetometer = magnetometer / norm(magnetometer);
            
            %quaternion ingegration
            half_dt = -0.5 * dt;
            w = [0; wx; wy; wz];
            q_dot = obj.quaternion_mult(w, obj.q_last);
            q_gyro = [obj.q_last(1) + q_dot(1) * half_dt;
                      obj.q_last(2) + q_dot(2) * half_dt;
                      obj.q_last(3) + q_dot(3) * half_dt;
                      obj.q_last(4) + q_dot(4) * half_dt];
            q_gyro = obj.quat_normalize(q_gyro);
            
            %calculate predicted gravity vector
            conj_q_gyro = obj.quaternion_conj(q_gyro);
            R_gyro = obj.prepare_body_to_earth_rotation_matrix(conj_q_gyro);
            g_predict = R_gyro * accelerometer;
            g_predict = g_predict / norm(g_predict);
                        
            %calculate predicted magnetic field vector
            l_predict = R_gyro * magnetometer; 
            l_predict = l_predict / norm(l_predict);
            
            %calculate delta change of quaternion for fusing gyroscope and aceelerometer
            weight_accel = 0.005;
            delta_q_acc = obj.convert_gravity_to_quat(g_predict);
            delta_q_acc = obj.quat_normalize(delta_q_acc);
            bar_delta_q_acc = obj.lerp(q_identy, delta_q_acc, weight_accel);
            bar_delta_q_acc = obj.quat_normalize(bar_delta_q_acc);
            
            %calculate delta change of quaternion for fusing gyroscope and magnetometer
            weight_mag = 0.005;
            delta_q_mag = obj.convert_magnetic_field_to_quat(l_predict);
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
            euler_angles = obj.quat_to_euler(obj.q_estimate);

            obj.R = obj.quat_to_rotation_matrix(obj.q_estimate);
            
            obj.roll = euler_angles(1);
            obj.pitch = euler_angles(2);
            obj.yaw = euler_angles(3);
            
            ret_obj = obj;
        end
    end
end