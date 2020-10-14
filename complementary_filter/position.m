classdef position
    properties
        home_longitude = 0;
        home_latitude = 0;
        home_ecef_x = 0;
        home_ecef_y = 0;
        home_ecef_z = 0;
        
        fused_enu_x = 0;
        fused_enu_y = 0;
        fused_enu_z = 0;
        fused_enu_vx = 0;
        fused_enu_vy = 0;
        fused_enu_vz = 0;
    end
    
    methods
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
        
        function ret_obj = filter_state_init(obj, longitude, latitude, height_ams, ...
                                             gps_ned_vx, gps_ned_vy, height_rate)
            enu_pos = obj.convert_gps_ellipsoid_coordinates_to_enu(longitude, latitude, height_ams);
            obj.fused_enu_x = enu_pos(1);
            obj.fused_enu_y = enu_pos(2);
            obj.fused_enu_z = height_ams;
            
            obj.fused_enu_vx = gps_ned_vy;
            obj.fused_enu_vy = gps_ned_vx;
            obj.fused_enu_vz = height_rate;
            ret_obj = obj;
        end
        
        function ret_obj = state_predict(obj, R_body_to_earth, accel_b_x, accel_b_y, accel_b_z, dt)            
            accel_i = R_body_to_earth * [accel_b_x; accel_b_y; accel_b_z];
            
            accel_i_x = -accel_i(2);
            accel_i_y = -accel_i(1);
            accel_i_z = -(accel_i(3) + 9.8);
            
            obj.fused_enu_vx = obj.fused_enu_vx + (accel_i_x * dt);
            obj.fused_enu_vy = obj.fused_enu_vy + (accel_i_y * dt);
            obj.fused_enu_vz = obj.fused_enu_vz + (accel_i_z * dt);
            
            obj.fused_enu_x = obj.fused_enu_x + (obj.fused_enu_vx * dt) + ...
                              (0.5 * accel_i_x * dt *dt);
            obj.fused_enu_y = obj.fused_enu_y + (obj.fused_enu_vy * dt) + ...
                              (0.5 * accel_i_y * dt *dt);
            obj.fused_enu_z = obj.fused_enu_z + (obj.fused_enu_vz * dt) + ...
                              (0.5 * accel_i_z * dt *dt);
                          
            ret_obj = obj;
        end
        
        function ret_obj = gps_update(obj, long_ref, lat_ref, vx_ned_ref, vy_ned_ref)
            velocity_enu_x = vy_ned_ref;
            velocity_enu_y = vx_ned_ref;
            
            position_enu = obj.convert_gps_ellipsoid_coordinates_to_enu(long_ref, lat_ref, 0);
            position_enu_x = position_enu(1);
            position_enu_y = position_enu(2);
            
            %velocity complementary filter
            weight_gps_vel = 0.4;
            obj.fused_enu_vx = (1 - weight_gps_vel) * obj.fused_enu_vx + weight_gps_vel * velocity_enu_x;
            obj.fused_enu_vy = (1 - weight_gps_vel) * obj.fused_enu_vx + weight_gps_vel * velocity_enu_y;
            
            %position complementary filter
            weight_gps_pos = 0.995;
            obj.fused_enu_x = (1 - weight_gps_pos) * obj.fused_enu_x + weight_gps_pos * position_enu_x;
            obj.fused_enu_y = (1 - weight_gps_pos) * obj.fused_enu_x + weight_gps_pos * position_enu_y;
            
            ret_obj = obj;
        end
        
        function ret_obj = barometer_update(obj, height_ref, vz_ref)
            weight_barometer_vel = 0.005;
            obj.fused_enu_vz = (1 - weight_barometer_vel) * obj.fused_enu_vz + weight_barometer_vel * vz_ref;
            
            weight_barometer_height = 0.005;
            obj.fused_enu_z = (1 - weight_barometer_height) * obj.fused_enu_z + weight_barometer_height * height_ref;
            
            ret_obj = obj;
        end
    end
end