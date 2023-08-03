% Adaptive Control - Simulation 3
% Masoud Pourghavam
% Student Number: 810601044
% Question 1-2 with color noise

%% --------------------------------------------- %%

% Clear workspace, close figures, and clear command window
clear all;
close all;
clc;

%%

format long
s = tf('s');
z = tf('z');
Ts = 1;
Gz = (z - 0.32) / ((z - 0.12) * (z - 0.47));

[numD, denD] = tfdata(Gz);
numD = cell2mat(numD);
denD = cell2mat(denD);

pole = pole(Gz);
zero = zero(Gz);

%% desired poles
Mp = 0.1;
zeta = ((log(Mp)^2) / (pi^2 + log(Mp)^2))^0.5;
ts = 2;
sigma = 4 / ts;
wn = sigma / zeta;
s1 = -zeta * wn + i * (wn * (1 - zeta^2)^0.5);
s2 = -zeta * wn - i * (wn * (1 - zeta^2)^0.5);

Pz1 = exp(s1 * Ts);
Pz2 = exp(s2 * Ts);

n = 4;
a1 = denD(2);
a0 = denD(3);

b1 = numD(2);
b0 = numD(3);

c1 = -0.5;
c0 = 0.06;

r1 = 1;
r0 = -0.85;
s1 = 1.047;
s0 = -0.12;
teta0 = [r1; r0; s1; s0; c1; c0];
N = 300;
sample_number = zeros(N, 1);
for i = 1:N
    sample_number(i, 1) = i;
end

% Generate white noise zero mean and unit var
var_noise = 0.01;
white_noise = sqrt(var_noise) * randn(N, 1);
white_noise = white_noise - mean(white_noise);

% Square pulse reference input
uc_in = zeros(N, 1);
for i = 1:N / 4
    uc_in(i, 1) = 1;
end
for i = N / 4:N / 2
    uc_in(i, 1) = -1;
end
for i = N / 2:3 * N / 4
    uc_in(i, 1) = 1;
end
for i = 3 * N / 4:N
    uc_in(i, 1) = -1;
end

% Calculate output without noise
y_out = zeros(N, 1);
epsilon = zeros(N, 1);
% Initial condition
phi_out = zeros(1, n);
k = 1;
sse = 0;
I = 0;
NN = 6; % number of controller parameter
teta_hat = [0.1; 0.1; 0.3; 0.3; 0.1; 0.1];
p = 1000 * eye(NN);
phi_t = zeros(1, NN);

r_hat1 = zeros(N, 1);
r_hat0 = zeros(N, 1);

s_hat1 = zeros(N, 1);
s_hat0 = zeros(N, 1);

c_hat1 = zeros(N, 1);
c_hat0 = zeros(N, 1);

u_control = zeros(N, 1);

RR = zeros(1, 1);
SS = zeros(1, 2);

Am = (z - Pz1) * (z - Pz2);
Am_coef = cell2mat(tfdata(Am));
degAm = 2;
degBm = 1;
degB_plus = degBm;
degB_minus = degBm - degB_plus;
d0 = 1;

bm0 = evalfr(Am, 1);
Bm = bm0 * z^(degBm);

deg_A0 = degAm - degB_plus - 1;
A0 = z^(d0 - 1);
A0_star = A0 * (z^(1 - d0));
T = A0_star * evalfr(Am, 1);
T_coef = cell2mat(tfdata(T));
t1_hat = T_coef(1);

yf = zeros(N, 1);
uf = zeros(N, 1);
ef = zeros(N, 1);
yfs = zeros(1, 2);
ufs = zeros(1, 2);
efs = zeros(1, 2);

for k = 2:N-1
    
    for ii = 1:2
        if k - ii <= 0
            ud = 0;
            yd = 0;
            ed = 0;
        else
            ud = -uf(k - ii, 1);
            yd = -yf(k - ii, 1);
            ed = -ef(k - ii, 1);
        end
        ufs(1, ii) = ud;
        yfs(1, ii) = yd;
        efs(1, ii) = ed;
    end
    uf(k, 1) = u_control(k, 1) + ufs * [Am_coef(2); Am_coef(3)];
    yf(k, 1) = y_out(k, 1) + yfs * [Am_coef(2); Am_coef(3)];
    ef(k, 1) = white_noise(k, 1) + efs * [Am_coef(2); Am_coef(3)];
   
    for m = 0:1
        if k - m - d0 <= 0
            y_dummy = 0;
            u_dummy = 0;
            e_dummy = 0;
        else 
            y_dummy = yf(k - m - d0, 1);
            u_dummy = uf(k - m - d0, 1);
            e_dummy = ef(k - m - d0, 1);
        end
        phi_t(1, m + 1) = u_dummy;
        phi_t(1, m + 3) = y_dummy;
        phi_t(1, m + 5) = e_dummy;
    end
    p = p - ((p * (phi_t') * phi_t * p) / (1 + phi_t * p * (phi_t')));
    gain = p * (phi_t');
    
    y_out(k, 1) = phi_t * teta0 + ef(k + 1 - d0, 1);
    epsilon(k, 1) = y_out(k, 1) - phi_t * teta_hat;
    teta_hat = teta_hat + gain * epsilon(k, 1);
    
    r_hat1(k, 1) = teta_hat(1);
    r_hat0(k, 1) = teta_hat(2);   

    s_hat1(k, 1) = teta_hat(3);
    s_hat0(k, 1) = teta_hat(4);
    
    c_hat1(k, 1) = teta_hat(5);
    c_hat0(k, 1) = teta_hat(6);
    
    for ii = 0:0
        if k - ii <= 0
            uu = 0;
        else
            uu = -u_control(k - ii, 1);
        end
        RR(1, ii + 1) = uu;
    end
    for ii = 0:1
        if k - ii <= 0
            yy = 0;
        else 
            yy = y_out(k - ii + 1, 1);
        end
        SS(1, ii + 1) = yy;    
    end
    
    u_control(k + 1, 1) = (RR * (r_hat0(k, 1)) + t1_hat * uc_in(k, 1) - SS * [s_hat1(k, 1); s_hat0(k, 1)]) / r_hat1(k, 1);
   
    I = (abs(y_out(k, 1) - (phi_t * teta_hat)));
    sse = I;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
plot(sample_number, uc_in, "black", 'LineWidth', 1.5);
hold on;
plot(sample_number, y_out, "red", 'LineWidth', 1.5);
ylim([-2 2]);
xlabel('Iter');
legend('Uc', 'y');

figure;
plot(sample_number, u_control, "black", 'LineWidth', 1.5);
ylim([-50 50]);
xlabel('Iter');
ylabel('u control');

figure;
subplot(2, 1, 1);
plot(sample_number, r_hat1, "black", 'LineWidth', 1.5);
xlabel('Iter');
ylabel('rhat1');
subplot(2, 1, 2);
plot(sample_number, r_hat0, "black", 'LineWidth', 1.5);
xlabel('Iter');
ylabel('rhat0');

figure;
subplot(2, 1, 1);
plot(sample_number, s_hat1, "black", 'LineWidth', 1.5);
xlabel('Iter');
ylabel('shat1');
subplot(2, 1, 2);
plot(sample_number, s_hat0, "black", 'LineWidth', 1.5);
xlabel('Iter');
ylabel('shat0');

variance_y = var(y_out);
mean_y = mean(y_out);
variance_Ucontrol = var(u_control);
mean_Ucontrol = mean(u_control);

output_variance = sprintf('Output variance is %d', variance_y);
output_mean = sprintf('Output mean is %d', mean_y);
control_variance = sprintf('Control signal variance is %d', variance_Ucontrol);
control_mean = sprintf('Control signal mean is %d', mean_Ucontrol);

output_variance, output_mean, control_variance, control_mean
