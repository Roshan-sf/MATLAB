%% Control Gain Selection Example

clear
clc

%% Example 17.1 from the book
%
% t_r <= 30 seconds
% M_p <= 30%
% t_s <= 100 seconds

clear
clc
% We found s = -.05 +- j*0.1
kp = 0.0125;
kd = .1;
I = 1;
wn = sqrt(kp/I)
zeta = kd/(2*wn*I)

Gsc = tf(1, [I 0 0])
step(Gsc)

Gc_d = tf([kd 0], 1)
Gp = feedback(Gsc, Gc_d)
Gc = tf(kp, 1)

T = feedback(Gc*Gp, 1)

figure
step(T)

%% Single Axis Analysis
Mp_reg = .1;
ts = 100;
% tr < 12;

zeta = sqrt(log(Mp_reg)^2/(pi^2 + log(Mp_reg)^2));

wn = log(0.02*sqrt(1-zeta^2))/-zeta/ts;

beta = atan(sqrt(1-zeta^2)/zeta);
tr = (pi-beta)/wn/sqrt(1-zeta^2);

% Extend to each Axis

J = diag([20 30 40]);

Kp = 2*J*eye(3)*wn^2
Kd = J*eye(3)*2*zeta*wn

epsilon_b_ECI_0 = [-.2; .4; .2];
q_b_ECI_0 = [epsilon_b_ECI_0; sqrt(1-norm(epsilon_b_ECI_0)^2)];
w_0 = [.1; -.1; .2];

%epsilon_C = [.1; -.3; .4];
epsilon_C = [0; 0; 0];
q_C = [epsilon_C; sqrt(1-norm(epsilon_C)^2)];

tspan = 120;
out = sim('gain_selection_ex');

subplot(2,1,1)
plot(q_b_ECI.time, squeeze(q_b_ECI.signals.values))
xlabel('Time (seconds)')
ylabel('Quaternion Components')
legend('\epsilon_x', '\epsilon_y', '\epsilon_z', '\eta')
grid on

subplot(2,1,2)
plot(omega.time, squeeze(omega.signals.values))
xlabel('Time (seconds)')
ylabel('Body Rates (rad/sec)')
legend('\omega_x', '\omega_y', '\omega_z')
grid on