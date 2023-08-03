% Adaptive Control - Simulation 3
% Masoud Pourghavam
% Student Number: 810601044
% Question 1-1 with color noise

%% Initialization %%
% Clear workspace, close figures, and clear command window
clear all;
close all;
clc;

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

% Desired pole location
Mp = 0.1;
zeta = ((log(Mp)^2) / (pi^2 + log(Mp)^2))^0.5;
ts = 2;
sigma = 4 / ts;
wn = sigma / zeta;
s1 = -zeta * wn + i * (wn * (1 - zeta^2)^0.5);
s2 = -zeta * wn - i * (wn * (1 - zeta^2)^0.5);

Pz1 = exp(s1 * Ts);
Pz2 = exp(s2 * Ts);

Am = (z - Pz1) * (z - Pz2);

n = 4;
a1 = denD(2);
a0 = denD(3);

b1 = numD(2);
b0 = numD(3);

teta0 = [a1; a0; b1; b0];
N = 100;
sample_number = zeros(N, 1);
for i = 1:N
    sample_number(i, 1) = i;
end

%% Square pulse reference input
uc_in = zeros(N, 1);
for i = 1:N/4
    uc_in(i, 1) = 1;
end
for i = N/4:N/2
    uc_in(i, 1) = -1;
end
for i = N/2:3*N/4
    uc_in(i, 1) = 1;
end
for i = 3*N/4:N
    uc_in(i, 1) = -1;
end

%% Calculate output without noise
y_out = zeros(N, 1);

% Initial condition
a_hat1 = zeros(N, 1);
a_hat0 = zeros(N, 1);

b_hat1 = zeros(N, 1);
b_hat0 = zeros(N, 1);

c_hat1 = zeros(N, 1);
c_hat0 = zeros(N, 1);

r0 = zeros(N, 1);
t0 = zeros(N, 1);

s1 = zeros(N, 1);
s0 = zeros(N, 1);
u_control = zeros(N, 1);
RR = zeros(1, 1);
SS = zeros(1, 2);

epsilon = zeros(N, 1);

% Generate white noise
var_noise = 0.01;
white_noise = sqrt(var_noise) * randn(N, 1);
white_noise = white_noise - mean(white_noise);

c1 = -0.48;
c0 = 0.071;

teta0 = [teta0; c1; c0];

% Initial condition
NN = 6;
teta_hat = [0.5; 0.5; 0.5; 1; 1; 0.05];
p = 1000 * eye(NN);
phi_t = zeros(1, NN);
k = 1;
sse = 10;
I = 0;

degA_hat = 2;
degB_hat = 1;
degB_plus = degB_hat;
degB_minus = degB_hat - degB_plus;
degBm = degB_hat;
bm0 = evalfr(Am, 1);
Bm = bm0 * z^(degBm);
deg_Bprim_m = degBm - degB_minus;
deg_A0 = degA_hat - degB_plus - 1;
R_prim = 1;
A0 = z^(deg_A0);

for k = 1:N
    for m = 1:2
        if k - m <= 0
            y_dummy = 0;
            u_dummy = 0;
            e_dummy = 0;
        else
            y_dummy = -y_out(k - m, 1);
            u_dummy = u_control(k - m, 1);
            e_dummy = white_noise(k - m, 1);
        end
        phi_t(1, m) = y_dummy;
        phi_t(1, m + 2) = u_dummy;
        phi_t(1, m + 4) = e_dummy;
    end
    p = p - ((p * (phi_t') * phi_t * p) / (1 + phi_t * p * (phi_t')));
    gain = p * (phi_t');
    y_out(k, 1) = phi_t * teta0 + white_noise(k, 1);
    epsilon(k, 1) = y_out(k, 1) - phi_t * teta_hat;
    teta_hat = teta_hat + gain * (epsilon(k, 1));

    a_hat1(k, 1) = teta_hat(1);
    a_hat0(k, 1) = teta_hat(2);

    b_hat1(k, 1) = teta_hat(3);
    b_hat0(k, 1) = teta_hat(4);

    c_hat1(k, 1) = teta_hat(5);
    c_hat0(k, 1) = teta_hat(6);

    A_hat = z^2 + a_hat1(k, 1) * z + a_hat0(k, 1);
    B_hat = b_hat1(k, 1) * z + b_hat0(k, 1);

    B_minus = b_hat1(k, 1);
    B_plus = [b_hat1(k, 1) / B_minus, b_hat0(k, 1) / B_minus];

    Bprim_m = (bm0 / B_minus) * z^(deg_Bprim_m);

    R = R_prim * B_plus;

    r0(k, 1) = R(2);

    T = A0 * Bprim_m;
    Tcoef = cell2mat(tfdata(T));
    t0(k, 1) = Tcoef(1);
    S = minreal((A0 * Am - A_hat) / B_minus);
    Scoef = cell2mat(tfdata(S));

    s1(k, 1) = Scoef(1);
    s0(k, 1) = Scoef(2);

    for ii = 1:1
        if k - ii <= 0
            uu = 0;
        else
            uu = -u_control(k - ii, 1);
        end
        RR(1, ii) = uu;
    end
    for ii = 0:1
        if k - ii <= 0
            yy = 0;
        else
            yy = y_out(k - ii, 1);
        end
        SS(1, ii + 1) = yy;
    end

    u_control(k, 1) = RR * r0(k, 1) + t0(k, 1) * uc_in(k, 1) - SS * [s1(k, 1); s0(k, 1)];

    I = (abs(y_out(k, 1) - (phi_t * teta_hat)));
    sse = I;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure 
plot(sample_number,uc_in,"black",'LineWidth',1.5)
hold on 
plot(sample_number,y_out,"red",'LineWidth',1.5)
ylim([-2 2])
xlabel('Iter')
legend('Uc','y')
figure
plot(sample_number,u_control,"black",'LineWidth',1.5)
ylim([-40 40])
xlabel('Iter')
ylabel('u control')

%%
figure
subplot(2,1,1)
plot(sample_number,a_hat1,"black",'LineWidth',1.5)
xlabel('Iter')
ylabel('ahat1')
subplot(2,1,2)
plot(sample_number,a_hat0,"black",'LineWidth',1.5)
xlabel('Iter')
ylabel('ahat0')
%%
figure
subplot(2,1,1)
plot(sample_number,b_hat1,"black",'LineWidth',1.5)
xlabel('Iter')
ylabel('bhat1')
subplot(2,1,2)
plot(sample_number,b_hat0,"black",'LineWidth',1.5)
xlabel('Iter')
ylabel('bhat0')
%%
figure
subplot(2,1,1)
plot(sample_number,c_hat1,"black",'LineWidth',1.5)
xlabel('Iter')
ylabel('chat1')
subplot(2,1,2)
plot(sample_number,c_hat0,"black",'LineWidth',1.5)
xlabel('Iter')
ylabel('chat0')

%%
figure
subplot(1,1,1)
plot(sample_number,r0,"black", 'LineWidth',1.5)
xlabel('Iter')
ylabel('r0')

%%
figure
subplot(2,1,1)
plot(sample_number,s0,"black",'LineWidth',1.5)
xlabel('Iter')
ylabel('s0')
subplot(2,1,2)
plot(sample_number,s1,"black",'LineWidth',1.5)
xlabel('Iter')
ylabel('s1')

%%
figure
plot(sample_number,t0,"black",'LineWidth',1.5)
xlabel('Iter')
ylabel('t0')

%%
variance_y = var(y_out);
mean_y = mean(y_out);
variance_Ucontrol = var(u_control);
mean_Ucontrol = mean(u_control);

sprintf('output variance = %d', variance_y)
sprintf('output mean = %d', mean_y)
sprintf('control signal variance = %d', variance_Ucontrol)
sprintf('signal control mean = %d', mean_Ucontrol)
