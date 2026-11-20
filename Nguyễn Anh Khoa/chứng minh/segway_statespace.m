%% segway_statespace.m
% Chuyen phuong trinh vi phan (35) va (38) trong PDF sang khong gian trang thai
% Ky hieu theo PDF: M, m, Iw, Iy, Iz, R, l, D, g
close all; clear; clc;

%% 1. Thong so (lay tu code cu cua ban)
M  = 4;        % khoi luong than [kg]
m  = 0.3;      % khoi luong rotor banh [kg]
Iw = 0.000735; % quan tinh banh [kg m^2]
Iy = 0.0377;   % quan tinh than quanh truc y (pitch)
Iz = 0.0406;   % quan tinh than quanh truc z (yaw)
R  = 0.08;     % ban kinh banh [m]
D  = 0.29;     % khoang cach 2 banh [m]
l  = 0.33;     % khoang cach tam khoi luong - truc banh [m]
g  = 9.81;

%% 2. He phang (tien/lui + nghieng): x = [x; xdot; theta; thetadot], u = TL+TR
a     = 2*Iw/R^2 + 2*m + M;     % he so "khoi luong tuong duong" o pt (35) dong 1
c     = Iy + M*l^2;             % he so o pt (35) dong 2
Delta = a*c - (M*l)^2;          % dinh thuc ma tran khoi luong

% Giai xdd, thetadd tu pt (35):
%   a*xdd   + M*l*thdd = u/R
%   M*l*xdd + c*thdd   = M*g*l*theta - u
xdd_theta = -(M*l)*(M*g*l)/Delta;
xdd_u     = (c/R + M*l)/Delta;
thdd_theta =  a*(M*g*l)/Delta;
thdd_u     = -(M*l/R + a)/Delta;

Ap = [0 1 0 0;
      0 0 xdd_theta 0;
      0 0 0 1;
      0 0 thdd_theta 0];
Bp = [0; xdd_u; 0; thdd_u];
Cp = eye(4);
Dp = zeros(4,1);

sysP = ss(Ap, Bp, Cp, Dp);
sysP.StateName = {'x','xdot','theta','thetadot'};
sysP.InputName = {'TL+TR'};

%% 3. He quay (yaw): z = [psi; psidot], v = TL-TR  (pt 38)
Jpsi = R*(2*Iz/D + Iw*D/R^2 + m*D);
Ay = [0 1; 0 0];
By = [0; 1/Jpsi];
sysY = ss(Ay, By, eye(2), zeros(2,1));
sysY.StateName = {'psi','psidot'};
sysY.InputName = {'TL-TR'};

%% 4. Kiem tra he thong
disp('--- Ma tran A, B (he phang) ---'); disp(Ap); disp(Bp);
disp('--- Cuc cua he phang (se co 1 cuc duong = mat on dinh) ---');
disp(eig(Ap));
disp(['Hang ma tran dieu khien (can = 4): ', num2str(rank(ctrb(Ap,Bp)))]);

%% 5. Thiet ke LQR: u = -K x
Q = diag([1, 1, 100, 1]);   % phat sai lech x, xdot, theta, thetadot
Rw = 1;                      % phat nang luong dieu khien
K = lqr(Ap, Bp, Q, Rw);
disp('K = '); disp(K);

sysCL = ss(Ap - Bp*K, Bp, Cp, Dp);
disp('--- Cuc vong kin (phai nam ben trai) ---'); disp(eig(sysCL.A));

% Dap ung: robot bi nghieng 0.1 rad luc dau
x0 = [0; 0; 0.1; 0];
figure;
[y, t] = initial(sysCL, x0, 5);
subplot(2,1,1); plot(t, y(:,3)*180/pi, 'LineWidth', 2); grid on
ylabel('theta [do]'); title('Robot nghieng 0.1 rad luc dau, bo dieu khien keo ve');
subplot(2,1,2); plot(t, y(:,1), 'LineWidth', 2); grid on
ylabel('x [m]'); xlabel('t [s]');

%% 6. Quan sat anh huong cua l (ly do can gain scheduling)
ls = linspace(0.2598, 0.4281, 5);
figure; hold on
for i = 1:numel(ls)
    li = ls(i);
    ci = Iy + M*li^2; Di = a*ci - (M*li)^2;
    Ai = [0 1 0 0;
          0 0 -(M*li)*(M*g*li)/Di 0;
          0 0 0 1;
          0 0  a*(M*g*li)/Di 0];
    Bi = [0; (ci/R + M*li)/Di; 0; -(M*li/R + a)/Di];
    % dung K co dinh tinh o l=0.33 cho cac l khac
    p = eig(Ai - Bi*K);
    plot(real(p), imag(p), 'x', 'MarkerSize', 10, 'LineWidth', 2);
    leg{i} = sprintf('l = %.3f', li);
end
grid on; xlabel('Re'); ylabel('Im'); legend(leg);
title('Cuc vong kin khi dung K co dinh cho cac l khac nhau');


