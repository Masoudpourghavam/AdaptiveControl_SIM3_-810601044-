% Adaptive Control - Simulation 3
% Masoud Pourghavam
% Student Number: 810601044
% Question 1-4 Non-Adaptive MV

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

B=(z-0.32);
Bcoef=cell2mat(tfdata(B));

C=(z-0.6)*(z-0.17);
Ccoef=cell2mat(tfdata(C));

N = 500;
sample_number=zeros(N,1);

for i=1:N
    sample_number(i,1)=i;
end

y_out=zeros(N,1);
u_control=zeros(N,1);
V_out=zeros(N,1);
V_control=zeros(N,1);

% Generate white noise zero mean
var_noise=0.01;
white_noise=sqrt(var_noise)*randn(N,1);
white_noise=white_noise- mean(white_noise);

%%

n=2;
m=1;
d0=n-m;
F=1;
R=B;
S=C-A;
Scoef=cell2mat(tfdata(S));
s1=Scoef(1);
s0=Scoef(2);
Rcoef=cell2mat(tfdata(R));
r0=Rcoef(2);
sum=0;
yy=zeros(2,1);
uu=zeros(2,1);
ee=zeros(3,1);
%%
for k=1:N-1
    for ii=1:2
        if k-ii<=0
            y_dummy=0;
            u_dummy=0;
        else
            u_dummy=u_control(k-ii,1);
            y_dummy=-y_out(k-ii,1);
        end
        yy(ii,1)=y_dummy;
        uu(ii,1)=u_dummy;
    end
    for jj=0:2
        if k-jj<=0
            e_dummy=0;
        else
            e_dummy=white_noise(k-jj,1);
        end
        ee(jj+1,1)=e_dummy;
    end
    y_out(k,1)=[Acoef(2) Acoef(3)]*yy+Bcoef*uu+Ccoef*ee;
    if k <=1
        u_control(k,1)=0;
    elseif k>=2
        u_control(k,1)=-r0*u_control(k-1,1)-s1*y_out(k,1)-s0*y_out(k-1,1);
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
legend('y_out')
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