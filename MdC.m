clc;clear all;close all;
Tr=295;%环境温度
L=0.6;
Asf=0.2; 	%单位长度固体流体接触面积
Ac=0.04;%通道截面积
lambdas=10.6/10;
lambdaf=0.6/10;	
cps=Tr;
cpf=4180;	
rhof=1000;
rhos=7901;
muf=1.04e-3;
vx=0.1;
Pr=10;
epl=0.32;
d=5E-4;
Dh=2*d*epl/3/(1-epl);
Re=rhof*Dh*vx/muf;
hsf=(1-epl)/epl*lambdaf/Dh*(0.5*Re^0.5+0.2*Re^(2/3))*Pr^(1/3);	
f=0.16;
BT=1.5;
deltaTad=@(Ts)5.8-(295-Ts)/(293-280)*(5.8-3.9);
Tp=1.4;%流动时间，冷到热
Tq=1.4;%流动时间，热到冷
Tv=(1/f-Tp-Tq)/2;%时间，加磁
Tw=(1/f-Tp-Tq)/2;%时间，退磁
Numr=500;
%% 微分方程
% 数据容器与求解器
n=10;deltax=L/n;x=0:L/n:L;
m=10;deltat=Tp/m;
U=zeros(m+1,(n+1)*2);% 临时容器，存储半个循环的数据
Tf=zeros(2*(Numr+1)+1,n+1);
Ts=zeros(2*(Numr+1)+1,n+1);
theta=0.5;
% 定解条件
Tf0=Tr;Ts0=Tr;
% 初始化
Tf(1,:)=Tf0;Ts(1,:)=Ts0;
% 系数
A=epl*rhof*cpf;
B=epl*rhof*cpf*vx;
C=-lambdaf;
D=-hsf*Asf/Ac;
E=(1-epl)*rhos*cps;
F=-lambdas;
G=-hsf*Asf/Ac;
% 矩阵形成
Kp=zeros(2*(n+1));p=zeros(2*(n+1),1);
Kq=zeros(2*(n+1));q=zeros(2*(n+1),1);
% Tf Ts Tf Ts Tf Ts Tf Ts 
% 冷到热
Kp(1,1)=1;Kp(1,3)=-1;
Kp(2,2)=1;Kp(2,4)=-1;
for i=3:2:2*n
  Kp(i,i-2)=theta*C/deltax^2-theta*B/deltax;
  Kp(i,i)=A/deltat+theta*B/deltax-2*theta*C/deltax^2-D;
  Kp(i,i+2)=theta*C/deltax^2;
  Kp(i,i+1)=D;
  Kp(i+1,i-1)=theta*F/deltax^2;
  Kp(i+1,i+1)=E/deltat-2*theta*F/deltax^2-G;
  Kp(i+1,i+3)=theta*F/deltax^2;
  Kp(i+1,i)=G;
end
Kp(2*n+1,2*n-1)=1;Kp(2*n+1,2*n+1)=-1;
Kp(2*n+2,2*n)=1;Kp(2*n+2,2*n+2)=-1;
% 热到冷
B=-epl*rhof*cpf*vx;
deltat=Tq/m;
Kq(1,1)=1;Kq(1,3)=-1;
Kq(2,2)=1;Kq(2,4)=-1;
for i=3:2:2*n
  Kq(i,i-2)=theta*C/deltax^2;
  Kq(i,i)=A/deltat-theta*B/deltax-2*theta*C/deltax^2-D;
  Kq(i,i+2)=theta*B/deltax+theta*C/deltax^2;
  Kq(i,i+1)=D;
  Kq(i+1,i-1)=theta*F/deltax^2;
  Kq(i+1,i+1)=E/deltat-2*theta*F/deltax^2-G;
  Kq(i+1,i+3)=theta*F/deltax^2;
  Kq(i+1,i)=G;
end
Kq(2*n+1,2*n+1)=1;
Kq(2*n+2,2*n)=1;Kq(2*n+2,2*n+2)=-1;
%% 时间推进
counter=0;
while counter<=Numr
  % 加磁升温
  U(1,1:2:(2*n+1))=Tf(2*counter+1,:);
  U(1,2:2:(2*n+2))=Ts(2*counter+1,:)+deltaTad(Ts(2*counter+1,:));
  % 冷到热
  deltat=Tp/m;B=epl*rhof*cpf*vx;
  for j=2:(m+1)
    p(1)=0;p(2)=0;p(2*n+1)=0;p(2*n+2)=0;
    p(3:2:(2*n-1))=A*U(j-1,3:2:(2*n-1))/deltat-...
    (1-theta)*(B*(U(j-1,5:2:(2*n+1))-U(j-1,3:2:(2*n-1)))/deltax+...
    C*(U(j-1,5:2:(2*n+1))-2*U(j-1,3:2:(2*n-1))+U(j-1,1:2:(2*n-3)))/deltax^2);
    p(4:2:(2*n))=E*U(j-1,4:2:(2*n))/deltat-...
    (1-theta)*F*...
    (U(j-1,6:2:(2*n+2))-2*U(j-1,4:2:(2*n))+U(j-1,2:2:(2*n-2)))/deltax^2;
    U(j,:)=(Kp\p);
  end
  Tf(2*counter+2,:)=U(m+1,1:2:(2*n+1));
  Ts(2*counter+2,:)=U(m+1,2:2:(2*n+2));
  % 退磁降温
  U(1,1:2:(2*n+1))=U(m+1,1:2:(2*n+1));
  U(1,2:2:(2*n+2))=U(m+1,2:2:(2*n+2))-deltaTad(U(m+1,2:2:(2*n+2)));
  % 热到冷
  deltat=Tq/m;B=-epl*rhof*cpf*vx;
  for j=2:(m+1)
    q(1)=0;q(2)=0;q(2*n+1)=Tr;q(2*n+2)=0;
    q(3:2:(2*n-1))=A*U(j-1,3:2:(2*n-1))/deltat-...
    (1-theta)*(B*(U(j-1,5:2:(2*n+1))-U(j-1,3:2:(2*n-1)))/deltax+...
    C*(U(j-1,5:2:(2*n+1))-2*U(j-1,3:2:(2*n-1))+U(j-1,1:2:(2*n-3)))/deltax^2);
    q(4:2:(2*n))=E*U(j-1,4:2:(2*n))/deltat-...
    (1-theta)*F*...
    (U(j-1,6:2:(2*n+2))-2*U(j-1,4:2:(2*n))+U(j-1,2:2:(2*n-2)))/deltax^2;
    U(j,:)=(Kq\q);
  end
  Tf(2*counter+3,:)=U(m+1,1:2:(2*n+1));
  Ts(2*counter+3 ,:)=U(m+1,2:2:(2*n+2));
  counter=counter+1;
  counter
end
timeset=1/f/2*(1:length(Tf));
plot(timeset,Tf(:,1))


