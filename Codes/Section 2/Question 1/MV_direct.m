% Adaptive Control - Simulation 3
% Masoud Pourghavam
% Student Number: 810601044
% Question 2-1 MV Direct

%% --------------------------------------------- %%

% Clear workspace, close figures, and clear command window
clear all;
close all;
clc;

%%


format long



s=tf('s');
z=tf('z');

A=(z^2-1.2*z+0.5);
Acoef=cell2mat(tfdata(A));

B=(z-0.85);
Bcoef=cell2mat(tfdata(B));

C=z^2+0.5*z-0.25;
Ccoef=cell2mat(tfdata(C));

N = 1000;
sample_number=zeros(N,1);

for i=1:N
    sample_number(i,1)=i;
end

y_out=zeros(N,1);
epsilon=zeros(N,1);
u_control=zeros(N,1);
V_out=zeros(N,1);
V_control=zeros(N,1);

% Generate white noise zero mean
var_noise=0.01;
white_noise=sqrt(var_noise)*randn(N,1);
white_noise=white_noise- mean(white_noise);

%

a1=Acoef(2);
a0=Acoef(3);
b1=Bcoef(1);
b0=Bcoef(2);
c1=Ccoef(2);
c0=Ccoef(3);

nn=6;   % number of system parameter
teta0=[a1;a0;b1;b0;c1;c0];
NN=4;
teta_hat=[0.1;0.1;0.1;0.1];
p=(0.001)*eye(NN);
phi_t=zeros(1,NN);
sse=0.10;
I=0;

r_hat1=zeros(N,1);
r_hat0=zeros(N,1);
s_hat1=zeros(N,1);
s_hat0=zeros(N,1);
n=2;
m=1;
d0=n-m;
F=1;
d=d0;

%

yf=zeros(N,1);
uf=zeros(N,1);
phi_u=zeros(1,2);
phi_y=zeros(1,2);
phi_out=zeros(1,nn);

duty_cycle = 0.8;  % Increase the duty cycle (range: 0 to 1)
d = round(duty_cycle * d0);

for k=2:N-1
    
    for jj=1:2
        if k-jj<=0
            y_out_dummy=0;
            u_out_dummy=0;
            e_dummy=0;
        else
            u_out_dummy=u_control(k-jj,1);
            y_out_dummy=-y_out(k-jj,1);
            e_dummy=white_noise(k-jj,1);
        end
        phi_out(1,jj)=y_out_dummy;
        phi_out(1,jj+2)=u_out_dummy;
        phi_out(1,jj+4)=e_dummy;
    end
    
    y_out(k,1)=phi_out*teta0+white_noise(k,1);
    
    for m=1:2
        if k-m<=0
            uf_dummy=0;
            yf_dummy=0;
        else
            uf_dummy=-uf(k-m,1);
            yf_dummy=-yf(k-m,1);
        end
        phi_u(1,m)=uf_dummy;
        phi_y(1,m)=yf_dummy;
    end
    
    uf(k,1)=u_control(k,1)+phi_u*[c1;c0];
    yf(k,1)=y_out(k,1)+phi_y*[c1;c0];
    
    for ii=0:1
        if k-ii-d<=0
            y_dummy=0;
            u_dummy=0;
        else
            u_dummy=uf(k-ii-d,1);
            y_dummy=yf(k-ii-d,1);
        end
        phi_t(1,ii+1)=u_dummy;
        phi_t(1,ii+3)=y_dummy;
    end
    
    p=p-((p*(phi_t')*phi_t*p)/(1+phi_t*p*(phi_t')));
    gain=p*(phi_t');
    epsilon(k,1)=y_out(k,1)-phi_t*teta_hat;
    teta_hat=teta_hat+gain*(epsilon(k,1));
    
    r_hat1(k,1)=teta_hat(1);
    r_hat0(k,1)=teta_hat(2);   
    
    s_hat1(k,1)=teta_hat(3);
    s_hat0(k,1)=teta_hat(4);
    
    u_control(k+1,1)=(-r_hat0(k,1)*u_control(k,1)-s_hat1(k,1)*y_out(k,1)-s_hat0(k,1)*y_out(k-1,1))/r_hat1(k,1);

    if k == 2
        V_out(2,1) = y_out(2,1)^2;
        V_control(2,1) = u_control(2,1)^2;
    elseif k > 2
        V_out(k,1) = V_out(k-1,1) + y_out(k,1)^2;
        V_control(k,1) = V_control(k-1,1) + u_control(k,1)^2;
    end

end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(sample_number,y_out,"black",'LineWidth',1.5)
ylim([-0.5 0.5])
xlabel('Iter')
legend('y out')
figure 
plot(sample_number,u_control,"red",'LineWidth',1.5)
ylim([-0.5 0.5])
xlabel('Iter')
legend('Uc')

figure
plot(sample_number,V_out,"black",'LineWidth',1.5)
xlabel('Iter')
legend('accumulated loss for y')

figure
plot(sample_number,V_control,"red",'LineWidth',1.5)
xlabel('Iter')
legend('accumulated loss for Uc')
%%
figure
subplot(2,1,1)
plot(sample_number,r_hat1,"black",'LineWidth',1.5)
xlabel('Iter')
ylabel('rhat1')
subplot(2,1,2)
plot(sample_number,r_hat0,"black",'LineWidth',1.5)
xlabel('Iter')
ylabel('rhat0')
%
figure
subplot(2,1,1)
plot(sample_number,s_hat1,"black",'LineWidth',1.5)
xlabel('Iter')
ylabel('shat1')
subplot(2,1,2)
plot(sample_number,s_hat0,"black",'LineWidth',1.5)
xlabel('Iter')
ylabel('shat0')
%%
variance_y = var(y_out);
variance_noise = var(white_noise);
variance_Ucontrol = var(u_control);
mean_y = mean(y_out);
mean_noise = mean(white_noise);
mean_Ucontrol = mean(u_control);

fprintf('Output variance is %f\n', variance_y);
fprintf('Noise variance is %f\n', variance_noise);
fprintf('Control signal variance is %f\n', variance_Ucontrol);
fprintf('Output mean is %f\n', mean_y);
fprintf('Noise mean is %f\n', mean_noise);
fprintf('Control signal mean is %f\n', mean_Ucontrol);
