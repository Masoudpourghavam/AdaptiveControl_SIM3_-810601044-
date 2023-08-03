% Adaptive Control - Simulation 3
% Masoud Pourghavam
% Student Number: 810601044
% Question 1-5 Adaptive MV

%% --------------------------------------------- %%

% Clear workspace, close figures, and clear command window
clear all;
close all;
clc;

%%

format long


s=tf('s');
z=tf('z');

A=((z-0.12)*(z-0.47));
Acoef=cell2mat(tfdata(A));

B=(z-(1/0.32));
Bcoef=cell2mat(tfdata(B));

C=(z-0.6)*(z-0.17);
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

teta0=[a1;a0;b1;b0;c1;c0];
NN=6;
teta_hat=[1;1;1;1;1;1];
p=(10^5)*eye(NN);
phi_t=zeros(1,NN);
sse=10;
I=0;

a_hat1=zeros(N,1);
a_hat0=zeros(N,1);
b_hat1=zeros(N,1);
b_hat0=zeros(N,1);
c_hat1=zeros(N,1);
c_hat0=zeros(N,1);
% r1=zeros(N,1);
r0=zeros(N,1);
s1=zeros(N,1);
s0=zeros(N,1);


n=2;
m=1;
d0=n-m;
F=1;


%%
for k=1:N-1
    
    for ii=1:2
        if k-ii<=0
            y_dummy=0;
            u_dummy=0;
            e_dummy=0;
        else
            u_dummy=u_control(k-ii,1);
            y_dummy=-y_out(k-ii,1);
            e_dummy=epsilon(k-ii,1);
        end
        phi_t(1,ii)=y_dummy;
        phi_t(1,ii+2)=u_dummy;
        phi_t(1,ii+4)=e_dummy;
    end
    
    %p=p-((p*(phi_t')*phi_t*p)/(1+phi_t*p*(phi_t')));
    p=inv(inv(p)+(phi_t')*phi_t);
    gain=p*(phi_t');
    y_out(k,1)=phi_t*teta0+white_noise(k,1);
    epsilon(k,1)=y_out(k,1)-phi_t*teta_hat;
    teta_hat=teta_hat+gain*(epsilon(k,1));
    
    a_hat1(k,1)=teta_hat(1);
    a_hat0(k,1)=teta_hat(2);   
    
    b_hat1(k,1)=teta_hat(3);
    b_hat0(k,1)=teta_hat(4);
    
    c_hat1(k,1)=teta_hat(5);
    c_hat0(k,1)=teta_hat(6);
    
 
   
   r0(k,1)=(c_hat1(k,1)+(1/b_hat0(k,1))-a_hat1(k,1)-(c_hat0(k,1)/b_hat0(k,1))-(c_hat1(k,1)/(b_hat0(k,1)^2))+(a_hat0(k,1)/b_hat0(k,1))+(c_hat0(k,1)/(b_hat0(k,1)^3)))/(1+(a_hat0(k,1)/(b_hat0(k,1)^2))-(a_hat1(k,1)/b_hat0(k,1)));

   s0(k,1)=(c_hat0(k,1)/(b_hat0(k,1)^2))-(a_hat0(k,1)*r0(k,1)/b_hat0(k,1));

   s1(k,1)=(c_hat0(k,1)/b_hat0(k,1))+(c_hat1(k,1)/(b_hat0(k,1)^2))-(a_hat0(k,1)/b_hat0(k,1))-(s0(k,1)/b_hat0(k,1))-(a_hat1(k,1)*r0(k,1)/b_hat0(k,1));

    
    if k <=1
        u_control(k,1)=0;
    elseif k>=2
        u_control(k,1)=(-r0(k,1)*u_control(k-1,1)-s1(k,1)*y_out(k,1)-s0(k,1)*y_out(k-1,1));
    end
    if k<=1
        V_out(k,1)=y_out(k,1);
        V_control(k,1)=u_control(k,1);
    elseif k>=2
        V_out(k,1)=V_out(k-1,1)+y_out(k,1)^2;
        V_control(k,1)=V_control(k-1,1)+u_control(k,1)^2;
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
plot(sample_number,a_hat1,"black",'LineWidth',1.5)
xlabel('Iter')
ylabel('ahat1')
subplot(2,1,2)
plot(sample_number,a_hat0,"black",'LineWidth',1.5)
xlabel('Iter')
ylabel('ahat0')
%
figure
subplot(2,1,1)
plot(sample_number,b_hat1,"black",'LineWidth',1.5)
xlabel('Iter')
ylabel('bhat1')
subplot(2,1,2)
plot(sample_number,b_hat0,"black",'LineWidth',1.5)
xlabel('Iter')
ylabel('bhat0')
%
figure
subplot(2,1,1)
plot(sample_number,c_hat1,"black",'LineWidth',1.5)
xlabel('Iter')
ylabel('chat1')
subplot(2,1,2)
plot(sample_number,c_hat0,"black",'LineWidth',1.5)
xlabel('Iter')
ylabel('chat0')
%
figure
subplot(1,1,1)
plot(sample_number,r0,"black",'LineWidth',1.5)
xlabel('Iter')
ylabel('r0')
%
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
variance_y=var(y_out);
variance_noise=var(white_noise);
variance_Ucontrol=var(u_control);
mean_y=mean(y_out);
mean_noise=mean(white_noise);
mean_Ucontrol=mean(u_control);

sprintf('output variance is %d',variance_y)
sprintf('noise variance is %d',variance_noise)
sprintf('control signal variance is %d',variance_Ucontrol)
sprintf('output mean is %d',mean_y)
sprintf('noise mean is %d',mean_noise)
sprintf('signal control mean is %d',mean_Ucontrol)