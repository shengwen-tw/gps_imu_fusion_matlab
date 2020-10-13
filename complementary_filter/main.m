csv = csvread("gps_imu_compass_barometer_log1.csv");

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

ahrs = attitude;
ins = position;

%set home position
home_longitude = 120.9971619;
home_latitude = 24.7861614;
ins = ins.set_home_longitude_latitude(home_longitude, home_latitude, 0);

%record datas
roll = zeros(1, data_num);
pitch = zeros(1, data_num);
yaw = zeros(1, data_num);
gps_enu_x = zeros(1, data_num);
gps_enu_y = zeros(1, data_num);
gps_enu_z = zeros(1, data_num);

for i = 1: data_num
    ahrs = ...
        ahrs.complementary_filter(-accel_lpf_x(i), -accel_lpf_y(i), -accel_lpf_z(i), ...
                                  gyro_raw_x(i), gyro_raw_y(i), gyro_raw_z(i), ...
                                  mag_raw_x(i), mag_raw_y(i), mag_raw_z(i), dt);
    
    position_enu = ins.convert_gps_ellipsoid_coordinates_to_enu(longitude(i), latitude(i), barometer_height(i));
    
    gps_enu_x(i) = position_enu(1);
    gps_enu_y(i) = position_enu(2);
    gps_enu_z(i) = position_enu(3);
    
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

%2d position trajectory plot of x-y plane
figure('Name', 'x-y position (enu frame)');
plot(gps_enu_x, -gps_enu_y);
title('position (enu frame)');
xlabel('x [m]');
ylabel('y [m]');

%3d visualization of position trajectory
figure('Name', 'x-y-z position (enu frame)');
plot3(gps_enu_x, gps_enu_y, barometer_height);
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