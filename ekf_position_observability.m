pkg load control

dt = 0.2 %gps update rate = 5Hz

%state transition matrix
A = zeros(6, 6)
A(1:3, 1:3) = eye(3);
A(1:3, 4:6) = dt * eye(3);
A(4:6, 4:6) = eye(3);
disp(A)

%observation matrix
C = [1, 0, 0, 0, 0, 0;
     0, 1, 0, 0, 0, 0
     0, 0, 1, 0, 0, 0];

rank(obsv(A, C))
