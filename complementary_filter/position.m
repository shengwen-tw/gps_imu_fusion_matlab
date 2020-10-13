classdef position
    properties
        home_longitude;
        home_latitude;
        home_ecef_x;
        home_ecef_y;
        home_ecef_z;
    end
        
    methods
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
    end
end