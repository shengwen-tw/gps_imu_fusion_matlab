format long g
clear all

%csv = csvread("../dataset/fengyuan_20210705.csv");
csv = csvread("../dataset/nycu_engineer_building_fifth_20210707.csv");

%ms
timestamp_ms = csv(:, 1);
%m/s^2
accel_lpf_x = csv(:, 2);
accel_lpf_y = csv(:, 3);
accel_lpf_z = csv(:, 4);
%rad/s
gyro_raw_x = csv(:, 5);
gyro_raw_y = csv(:, 6);
gyro_raw_z = csv(:, 7);
%uT
mag_raw_x = csv(:, 8);
mag_raw_y = csv(:, 9);
mag_raw_z = csv(:, 10);
%degree
longitude = csv(:, 11);
latitude = csv(:, 12);
longitude = longitude .* 1e-7;
latitude = latitude .* 1e-7;
%m
gps_height_msl = csv(:, 13);
%m/s
gps_ned_vx = csv(:, 14);
gps_ned_vy = csv(:, 15);
gps_ned_vz = csv(:, 16);
%m
barometer_height = csv(:, 17);
%m/s
barometer_vz = csv(:, 18);

timestamp_s = timestamp_ms .* 0.001;
[data_num, dummy] = size(timestamp_ms);

%fusion period
dt = 0.01; %100Hz, 0.01s

ahrs = attitude_estimator;
ins = position_estimator;

%set home position
home_longitude = longitude(1);
home_latitude = latitude(1);
ins = ins.set_home_longitude_latitude(home_longitude, home_latitude, 0);

%record datas
roll = zeros(1, data_num);
pitch = zeros(1, data_num);
yaw = zeros(1, data_num);
%
gps_enu_x = zeros(1, data_num);
gps_enu_y = zeros(1, data_num);
gps_enu_z = zeros(1, data_num);
%
fused_enu_x = zeros(1, data_num);
fused_enu_y = zeros(1, data_num);
fused_enu_z = zeros(1, data_num);
%
fused_enu_vx = zeros(1, data_num);
fused_enu_vy = zeros(1, data_num);
fused_eny_vz = zeros(1, data_num);
%
%visualization of b1, b2, b3 vectors with quiver3
b1_visual_sample_cnt = 500;
quiver_cnt = floor(data_num / b1_visual_sample_cnt);
quiver_orig_x = zeros(1, quiver_cnt);
quiver_orig_y = zeros(1, quiver_cnt);
quiver_orig_z = zeros(1, quiver_cnt);
quiver_b1_u = zeros(1, quiver_cnt);
quiver_b1_v = zeros(1, quiver_cnt);
quiver_b1_w = zeros(1, quiver_cnt);
quiver_b2_u = zeros(1, quiver_cnt);
quiver_b2_v = zeros(1, quiver_cnt);
quiver_b2_w = zeros(1, quiver_cnt);
quiver_b3_u = zeros(1, quiver_cnt);
quiver_b3_v = zeros(1, quiver_cnt);
quiver_b3_w = zeros(1, quiver_cnt);
j = 1;
%corrected gravity vector
gravity_x_arr = zeros(1, data_num);
gravity_y_arr = zeros(1, data_num);
gravity_z_arr = zeros(1, data_num);
gravity_x_arr(1) = accel_lpf_x(1);
gravity_y_arr(1) = accel_lpf_y(2);
gravity_z_arr(1) = -accel_lpf_z(3);
%accelerometer norm and corrected norm
accelerometer_norm_arr = zeros(1, data_num);
gravity_norm_arr = zeros(1, data_num);

%ins state initialization
ins = ins.filter_state_init(longitude(1), latitude(1), barometer_height(1), ...
                            gps_ned_vx(1), gps_ned_vy(1), barometer_vz(1));

vel_ned_body = [0; 0; 0];
                        
for i = 2: data_num
    dt = timestamp_s(i) - timestamp_s(i - 1);
    
    %eliminate body-fixed frame acceleration induced by rotation for
    %getting better gravity vector
    accel_translational = cross([gyro_raw_x(i); gyro_raw_y(i); gyro_raw_z(i)], vel_ned_body);
    gravity = [accel_translational(1) - accel_lpf_x(i);
               accel_translational(2) - accel_lpf_y(i);
               accel_translational(3) - accel_lpf_z(i)];
    
    gravity_x_arr(i) = gravity(1);
    gravity_y_arr(i) = gravity(2);
    gravity_z_arr(i) = gravity(3);
    
    %attitude estimation
    ahrs = ...
        ahrs.complementary_filter(gravity(1), gravity(2), gravity(3), ...
                                  gyro_raw_x(i), gyro_raw_y(i), gyro_raw_z(i), ...
                                  mag_raw_x(i), mag_raw_y(i), mag_raw_z(i), dt);
    
    %position and velocity estimation
    ins = ins.state_predict(ahrs.R, ...
                            accel_lpf_x(i), accel_lpf_y(i), accel_lpf_z(i), ...
                            dt);

    %update z directional state from barometer
    ins = ins.barometer_update(barometer_height(i), barometer_vz(i));
    
    %update x-y directional state from gps if new gps data is obtained
    if (abs(gps_ned_vx(i) - gps_ned_vx(i - 1)) > 1e-2 || abs(gps_ned_vy(i) - gps_ned_vy(i - 1)) > 1e-2)
        ins = ins.gps_update(longitude(i), latitude(i), gps_ned_vx(i), gps_ned_vy(i));
    end
    
    vel_ned_body = ahrs.R * [ins.fused_enu_vy; ins.fused_enu_vx; -ins.fused_enu_vz];
    
    fused_enu_x(i) =  ins.fused_enu_x;
    fused_enu_y(i) = ins.fused_enu_y;
    fused_enu_z(i) = ins.fused_enu_z;
    fused_enu_vx(i) = ins.fused_enu_vx;
    fused_enu_vy(i) = ins.fused_enu_vy;
    fused_enu_vz(i) = ins.fused_enu_vz;

    %gps, barometer raw signal and transform it into enu frame for
    %visualization
    position_enu = ins.convert_gps_ellipsoid_coordinates_to_enu(longitude(i), latitude(i), barometer_height(i));
    gps_enu_x(i) = position_enu(1);
    gps_enu_y(i) = position_enu(2);
    gps_enu_z(i) = position_enu(3);
    
    %collect rotation vectors for visualization
    if mod(i, b1_visual_sample_cnt) == 0
        quiver_orig_x(j) = position_enu(1);
        quiver_orig_y(j) = position_enu(2);
        quiver_orig_z(j) = position_enu(3);
        quiver_b1_u(j) = ahrs.R(2, 1);
        quiver_b1_v(j) = ahrs.R(1, 1);
        quiver_b1_w(j) = -ahrs.R(3, 1);
        quiver_b2_u(j) = ahrs.R(2, 2);
        quiver_b2_v(j) = ahrs.R(1, 2);
        quiver_b2_w(j) = -ahrs.R(3, 2);
        quiver_b3_u(j) = ahrs.R(2, 3);
        quiver_b3_v(j) = ahrs.R(1, 3);
        quiver_b3_w(j) = -ahrs.R(3, 3);
        j = j + 1;
    end
    
    accelerometer_norm_arr(i) = sqrt(accel_lpf_x(i) * accel_lpf_x(i) + ...
                                     accel_lpf_y(i) * accel_lpf_y(i) + ...
                                     accel_lpf_z(i) * accel_lpf_z(i));
    gravity_norm_arr(i) = sqrt(gravity(1) * gravity(1) + ...
                               gravity(2) * gravity(2) + ...
                               gravity(3) * gravity(3));
    
    roll(i) = ahrs.roll;
    pitch(i) =ahrs.pitch;
    yaw(i) = ahrs.yaw;                   
end

%%%%%%%%
% Plot %
%%%%%%%%

%accelerometer
figure('Name', 'accelerometer');
subplot (3, 1, 1);
plot(timestamp_s, accel_lpf_x);
title('accelerometer');
xlabel('time [s]');
ylabel('ax [m/s^2]');
subplot (3, 1, 2);
plot(timestamp_s, accel_lpf_y);
xlabel('time [s]');
ylabel('ay [m/s^2]');
subplot (3, 1, 3);
plot(timestamp_s, accel_lpf_z);
xlabel('time [s]');
ylabel('az [m/s^2]');

%gyroscope
figure('Name', 'gyroscope');
subplot (3, 1, 1);
plot(timestamp_s, rad2deg(gyro_raw_x));
title('gyroscope');
xlabel('time [s]');
ylabel('wx [rad/s]');
subplot (3, 1, 2);
plot(timestamp_s, rad2deg(gyro_raw_y));
xlabel('time [s]');
ylabel('wy [rad/s]');
subplot (3, 1, 3);
plot(timestamp_s, rad2deg(gyro_raw_z));
xlabel('time [s]');
ylabel('wz [rad/s]');

%magnatometer
figure('Name', 'magnetometer');
subplot (3, 1, 1);
plot(timestamp_s, mag_raw_x);
title('magnetometer');
xlabel('time [s]');
ylabel('mx [uT]');
subplot (3, 1, 2);
plot(timestamp_s, mag_raw_z);
xlabel('time [s]');
ylabel('my [uT]');
subplot (3, 1, 3);
plot(timestamp_s, mag_raw_z);
xlabel('time [s]');
ylabel('mx [uT]');

%gps position
figure('Name', 'gps position (wgs84)');
subplot (3, 1, 1);
plot(timestamp_s, longitude);
title('gps position (wgs84)');
xlabel('time [s]');
ylabel('longitude [deg]');
subplot (3, 1, 2);
plot(timestamp_s, latitude);
xlabel('time [s]');
ylabel('latitude [deg]');
subplot (3, 1, 3);
plot(timestamp_s, gps_height_msl);
xlabel('time [s]');
ylabel('height MSL [m]');

%gps velocity
figure('Name', 'gps velocity');
subplot (3, 1, 1);
plot(timestamp_s, gps_ned_vx);
title('gps velocity');
xlabel('time [s]');
ylabel('vx [m/s]');
subplot (3, 1, 2);
plot(timestamp_s, gps_ned_vy);
xlabel('time [s]');
ylabel('my [m/s]');
subplot (3, 1, 3);
plot(timestamp_s, gps_ned_vz);
xlabel('time [s]');
ylabel('vz [m/s]');

%barometer
figure('Name', 'barometer');
subplot (2, 1, 1);
plot(timestamp_s, barometer_height);
title('barometer');
xlabel('time [s]');
ylabel('z [m]');
subplot (2, 1, 2);
plot(timestamp_s, barometer_vz);
xlabel('time [s]');
ylabel('vz [m/s]');

%estimated roll, pitch and yaw angle
figure('Name', 'euler angles');
subplot (3, 1, 1);
plot(timestamp_s, roll);
title('euler angles');
xlabel('time [s]');
ylabel('roll [deg]');
subplot (3, 1, 2);
plot(timestamp_s, pitch);
xlabel('time [s]');
ylabel('pitch [deg]');
subplot (3, 1, 3);
plot(timestamp_s, yaw);
xlabel('time [s]');
ylabel('yaw [deg]');

%position in enu frame
figure('Name', 'position (enu frame)');
subplot (3, 1, 1);
plot(timestamp_s, gps_enu_x);
title('position (enu frame)');
xlabel('time [s]');
ylabel('x [m]');
subplot (3, 1, 2);
plot(timestamp_s, gps_enu_y);
xlabel('time [s]');
ylabel('y [m]');
subplot (3, 1, 3);
plot(timestamp_s, barometer_height);
xlabel('time [s]');
ylabel('z [m]');

%fused velocity in enu frame
figure('Name', 'fused velocity (enu frame)');
grid on;
subplot (3, 1, 1);
hold on;
plot(timestamp_s, gps_ned_vy);
plot(timestamp_s, fused_enu_vx);
title('fused velocity (enu frame)');
legend('gps velocity', 'eskf velocity') ;
xlabel('time [s]');
ylabel('vx [m/s]');
subplot (3, 1, 2);
hold on;
plot(timestamp_s, gps_ned_vx);
plot(timestamp_s, fused_enu_vy);
xlabel('time [s]');
ylabel('vy [m/s]');
legend('gps velocity', 'eskf velocity') ;
subplot (3, 1, 3);
hold on;
plot(timestamp_s, barometer_vz);
plot(timestamp_s, fused_enu_vz);
legend('barometer velocity', 'eskf velocity') ;
xlabel('time [s]');
ylabel('vz [m/s]');

%accelerometer vs corrected gravity
figure('Name', 'gravity');
subplot (3, 1, 1);
hold on;
plot(timestamp_s, gravity_x_arr);
plot(timestamp_s, -accel_lpf_x);
title('gravity');
xlabel('time [s]');
ylabel('ax [m/s^2]');
legend('measured', 'corrected');
subplot (3, 1, 2);
hold on;
plot(timestamp_s, gravity_y_arr);
plot(timestamp_s, -accel_lpf_y);
xlabel('time [s]');
ylabel('ay [m/s^2]');
legend('measured', 'corrected');
subplot (3, 1, 3);
hold on;
plot(timestamp_s, gravity_z_arr);
plot(timestamp_s, -accel_lpf_z);
xlabel('time [s]');
ylabel('az [m/s^2]');
legend('measured', 'corrected');

%accelerometer norm vs corrected gravity norm
figure('Name', 'gravity norm');
hold on;
plot(timestamp_s(2: data_num), accelerometer_norm_arr(2: data_num));
plot(timestamp_s(2: data_num), gravity_norm_arr(2: data_num));
title('gravity norm');
xlabel('time [s]');
ylabel('ax (m/s^2]');
legend('measured', 'corrected');

%raw position vs fused position
figure('Name', 'raw position and fused position (enu frame)');
grid on;
subplot (3, 1, 1);
hold on;
plot(timestamp_s, gps_enu_x);
plot(timestamp_s, fused_enu_x);
title('raw position and fused position (enu frame)');
legend('gps x', 'eskf x') ;
xlabel('time [s]');
ylabel('x [m]');
subplot (3, 1, 2);
hold on;
plot(timestamp_s, gps_enu_y);
plot(timestamp_s, fused_enu_y);
xlabel('time [s]');
ylabel('y [m]');
legend('gps y', 'eskf y') ;
subplot (3, 1, 3);
hold on;
plot(timestamp_s, barometer_height);
plot(timestamp_s, fused_enu_z);
legend('barometer z', 'eskf z') ;
xlabel('time [s]');
ylabel('z [m]');

%2d position trajectory plot of x-y plane
figure('Name', 'x-y position (enu frame)');
grid on;
hold on;
plot(gps_enu_x, gps_enu_y, ...
     'Color', 'k', ...
     'Marker', 'o', ...
     'LineStyle', 'None', ...
     'MarkerSize', 3);
plot(fused_enu_x, fused_enu_y, 'Color', 'r');
legend('gps trajectory', 'fused trajectory') ;
title('position (enu frame)');
xlabel('x [m]');
ylabel('y [m]');

%3d visualization of position trajectory
figure('Name', 'x-y-z position (enu frame)');
hold on;
axis equal;
plot3(gps_enu_x, gps_enu_y, barometer_height, ...
      'Color', 'k', ...
      'Marker', 'o', ...
      'LineStyle', 'None', ...
      'MarkerSize', 2);
plot3(fused_enu_x, fused_enu_y, fused_enu_z, 'Color', 'm');
quiver3(quiver_orig_x,  quiver_orig_y,  quiver_orig_z, ...
        quiver_b1_u,  quiver_b1_v,  quiver_b1_w, 'Color', 'r', 'AutoScaleFactor', 0.2);
quiver3(quiver_orig_x,  quiver_orig_y,  quiver_orig_z, ...
        quiver_b2_u,  quiver_b2_v,  quiver_b2_w, 'Color', 'g', 'AutoScaleFactor', 0.2);
quiver3(quiver_orig_x,  quiver_orig_y,  quiver_orig_z, ...
        quiver_b3_u,  quiver_b3_v,  quiver_b3_w, 'Color', 'b', 'AutoScaleFactor', 0.2);
%legend('gps trajectory', 'b1 vector', 'b2 vector', 'b3 vector') 
legend('gps trajectory', 'fused trajectory', 'b1 vector', 'b2 vector', 'b3 vector') 
title('position (enu frame)');
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
daspect([1 1 1])
grid on
view(-10,20);

disp("Press any key to leave");
pause;
close all;