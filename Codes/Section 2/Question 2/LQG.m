% Adaptive Control - Simulation 3
% Masoud Pourghavam
% Student Number: 810601044
% Question 2-2 LQG

%% --------------------------------------------- %%

% Clear workspace, close figures, and clear command window
clear all;
close all;
clc;

%%


format long

s=tf('s');
z=tf('z');

B1=z-0.5;
B1coef=cell2mat(tfdata(B1));

b1_1=B1coef(1);
b1_0=B1coef(2);


A1=(z-0.3)*(z-0.45);
A1coef=cell2mat(tfdata(A1));
a1_1=A1coef(2);
a1_0=A1coef(3);


C1=z-0.7;
C1coef=cell2mat(tfdata(C1));
c1_1=C1coef(1);
c1_0=C1coef(2);


A2=z-0.9;
A2coef=cell2mat(tfdata(A2));
a2_1=A2coef(1);
a2_0=A2coef(2);



A=A1*A2;
Acoef=cell2mat(tfdata(A));
a2=Acoef(2);
a1=Acoef(3);
a0=Acoef(4);
n=3;



B=B1*A2;
Bcoef=cell2mat(tfdata(B));
b1=Bcoef(2);
b0=Bcoef(3);



C=A1*C1;
Ccoef=cell2mat(tfdata(C));
c2=Ccoef(2);
c1=Ccoef(3);
c0=Ccoef(4);



A2_plus=A2;
l=1;
A2_minus=1;
m=0;

%%

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

%square pulse refrence input
uc_in=zeros(N,1);


for i=1:N/4
    uc_in(i,1)=1;
end
for i=N/4:N/2
    uc_in(i,1)=-1;
end
for i=N/2:3*N/4
    uc_in(i,1)=1;
end
for i=3*N/4:N
    uc_in(i,1)=-1;
end


pho=2;
syms r p0 p1
eq1=r*p0-pho*a1_0;
eq2=r*(p1*p0+p1)-pho*(a1_1+a1_1*a1_0)-b1_0;
eq3=r*((p1^2)+(p0^2)+1)-pho*((a1_1^2)+(a1_0^2)+1)-(1+(b1_0^2));



eqs = [eq1, eq2, eq3];
[x,y,z]=vpasolve(eqs,[r,p1,p0]);

syms s
aa=double(vpa(solve(s^2+y(2)*s+z(2),s)));



p1=double(y(2));
p0=double(z(2));



r=double(x(2));

%%

syms r2 r1 r0 s3 s2 s1

eqq1=r2+a1_1+1-(c2+p1)==0;
eqq2=r1+a1_1*r2+a1_0+s2+1*b1_0-(c1+c2*p1+p0)==0;
eqq3=r0+a1_1*r1+a1_0*r2+s1+b1_0*s2-(c0+c1*p1+p0*c2)==0;
eqq4=a1_1*r0+a1_0*r1+b1_0*s1-(p1*c0+p0*c1)==0;
eqq5=a1_0*r0-(p0*c0)==0;
SOL=solve([eqq1,eqq2,eqq3,eqq4,eqq5],[r2,r1,r0,s2,s1]);

%%

r2=double(vpa(SOL.r2));
r1=double(vpa(SOL.r1));
r0=double(vpa(SOL.r0));
s2=double(vpa(SOL.s2));
s1=double(vpa(SOL.s1));

%%

T=((1+p1+p0)/(b1_1+b1_0))*C;
Tcoef=cell2mat(tfdata(T));
yy=zeros(3,1);
uu=zeros(3,1);
ee=zeros(4,1);
tt=zeros(4,1);
y_c=zeros(3,1);

%%

for k=2:N-1
    
    for ii=1:3
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
    for jj=0:3
        if k-jj<=0
            e_dummy=0;
        else
            e_dummy=white_noise(k-jj,1);
        end
        ee(jj+1,1)=e_dummy;
    end
    y_out(k,1)=[Acoef(2) Acoef(3) Acoef(4)]*yy+Bcoef*uu+Ccoef*ee;
   
    for ll=0:2
        if k-ll<=0
            yc_dummy=0;
        else 
            yc_dummy=y_out(k-ll,1);
        end
        y_c(ll+1,1)=yc_dummy;
    end
    for ll=0:3
        if k-ll<=0
            t_dummy=0;
        else 
            t_dummy=uc_in(k-ll,1);
        end
        tt(ll+1,1)=t_dummy;
    end
    
    u_control(k,1)=[-r2 -r1 -r0]*(uu)+[-1 -s2 -s1]*(y_c)+Tcoef*tt;
    
    if abs(u_control(k,1))>1000
         u_control(k,1)=1000*sign(u_control(k,1));
     end

    
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
%%
figure 
plot(sample_number,u_control,"black",'LineWidth',1.5)
ylim([-10 10])
xlabel('Iter')
ylabel('Control effort')
%%
mean_Ucontrol=mean(u_control);
mean_y=mean(y_out);
variance_Ucontrol=var(u_control);