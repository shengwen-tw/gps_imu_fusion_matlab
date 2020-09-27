pkg load symbolic

%===========================================%
% define math symbolic derivation variables %
%===========================================%

syms dt
syms ax_i ay_i az_i
syms px_last py_last pz_last vx_last vy_last vz_last
syms px_gps py_gps pz_gps vx_gps vy_gps vz_gps
syms p11 p22 p33 p44 p55 p66
syms q11 q22 q33 q44 q55 q66
syms r11 r22 r33 r44 r55 r66

x_last = [px_last; py_last; pz_last; vx_last; vy_last; vz_last];
u = [ax_i; ay_i; az_i];
z = [px_gps; py_gps; pz_gps; vx_gps; vy_gps; vz_gps];

%=====================================%
% define symbolic derivation matrices %
%=====================================%

P_last = [p11  0   0   0   0   0;
           0  p22  0   0   0   0;
           0   0  p33  0   0   0;
           0   0   0  p44  0   0;
           0   0   0   0  p55  0;
           0   0   0   0   0  p66];

Q = [q11   0   0   0   0   0;
       0 q22   0   0   0   0;
       0   0  q33  0   0   0;
       0   0   0  q44  0   0;
       0   0   0   0  q55  0;
       0   0   0   0   0  q66];

R = [r11   0   0   0   0   0;
       0 r22   0   0   0   0;
       0   0  r33  0   0   0;
       0   0   0  r44  0   0;
       0   0   0   0  r55  0;
       0   0   0   0   0  r66];

A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     0 0 0 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 0 0];

B = [0 0 0;
     0 0 0;
     0 0 0
     1 0 0
     0 1 0
     0 0 1];

H = [1 0 0 0 0 0;
     0 1 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 0 0];

I = eye(6);

%===========================%
% start symbolic derivation %
%===========================%

x_dot = A*x_last + B*u
x_priori = x_last + (dt .* x_dot)

P_dot = A*P_last + (P_last*transpose(A)) + Q
P_priori = P_last + (dt .* P_dot)

PHt = P_priori*H;
HPHt_R = H*P_priori*transpose(H) + R;
K = PHt*inv(HPHt_R)

x_posteriori = x_priori + K*(z - H*x_priori)
P_posteriori = (I - K*H)*P_priori
