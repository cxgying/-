clc;clear all;close all
%% 方程参数之定义
q=@(x,t)(0);%定义源项
a=0;b=1;L=b-a;%定义求解区域长度
time=0.6;%定义求解时间
%% 方程之精确解
U_exc=@(x,t)(cos(pi*t)*sin(pi*x));%精确解
%% 初边值条件之定义
%初值条件
g1=@(x)(sin(pi*x));
g2=@(x)(0);
dg1=@(x)(pi*cos(pi*x));
%狄利克雷边界
f1=@(t)(0);
f2=@(t)(0);
%% 配点
elementNum=100;%定义单元数
X=0:L/elementNum:L;
%% 求解器参数之定义
timestep=0.0001;%定义时间步长
xta=0.5;
timeNum=round(time/timestep);%计算时间步数
%% Solve
%% 基本系数之计算
alpha=0;beta=0;
k=timestep;h=(b-a)/elementNum;
a1=sin(h/2)*sin(h/2)/(sin(h)*sin(1.5*h));a2=2/(1+2*cos(h));a3=-3/(4*sin(1.5*h));
a4=3/(4*sin(1.5*h));a5=3*(1+3*cos(h))/(16*sin(0.5*h)*sin(0.5*h)*(2*cos(0.5*h)+cos(1.5*h)));
a6=-3*cos(0.5*h)*cos(0.5*h)/(2*sin(0.5*h)*sin(0.5*h)*(1+2*cos(h)));
W1=(1+2*alpha*k+k*k*xta*beta*beta)*a1-k*k*xta*a5;
W2=(1+2*alpha*k+k*k*xta*beta*beta)*a2-k*k*xta*a6;
W3=(2+2*alpha*k-(1-xta)*k*k*beta*beta)*a1+(1-xta)*k*k*a5;
W4=(2+2*alpha*k-(1-xta)*k*k*beta*beta)*a2+(1-xta)*k*k*a6;
% 待求解变量，每一行为一个时间步的数据
C=zeros(timeNum+1,elementNum+3);
%% 初始化
A=zeros(elementNum+3);
B=zeros(elementNum+3,1);
A(1,1)=a3;A(1,3)=a4;
B(1)=dg1(X(1));
for ictrl=2:elementNum+2
    A(ictrl,ictrl-1)=a1;
    A(ictrl,ictrl)=a2;
    A(ictrl,ictrl+1)=a1;
    B(ictrl)=g1(X(ictrl-1));
end
A(elementNum+3,elementNum+1)=a3;A(elementNum+3,elementNum+3)=a4;
B(elementNum+3)=dg1(X(elementNum+1));
C(1,:)=(A\B)';
%% 第一次时间步进
A=zeros(elementNum+3);
B=zeros(elementNum+3,1);
A(1,2)=a2+(-(W2+a2)*a1/(W1+a1));
B(1)=f1(k)+(-(W1+a1)*(W3*C(1,1)+W4*C(1,2)+W3*C(1,3)+2*k*g2(X(1))+k*k*q(X(1),k))/(W1+a1));
for ictrl=2:elementNum+2
    A(ictrl,ictrl-1)=W1+a1;
    A(ictrl,ictrl)=W2+a2;
    A(ictrl,ictrl+1)=W1+a1;
    B(ictrl)=W3*C(1,ictrl-1)+W4*C(1,ictrl)+W3*C(1,ictrl+1)+2*k*g2(X(ictrl-1))+k*k*q(X(ictrl-1),k);
end
A(elementNum+3,elementNum+2)=a2+(-(W2+a2)*a1/(W1+a1));
B(elementNum+3)=f2(k)+(-(W1+a1)*(W3*C(1,elementNum+1)+W4*C(1,elementNum+2)+W3*C(1,elementNum+3)+2*k*g2(X(elementNum+1))+k*k*q(X(elementNum+1),k))/(W1+a1));
C(2,:)=(A\B)';
%% 后续时间步进
A=zeros(elementNum+3);
A(1,2)=a2+(-a1*W2/W1);
for ictrl=2:elementNum+2
    A(ictrl,ictrl-1)=W1;
    A(ictrl,ictrl)=W2;
    A(ictrl,ictrl+1)=W1;
end
A(elementNum+3,elementNum+2)=a2+(-a1*W2/W1);
for jctrl=2:timeNum
    B(1)=f1(jctrl*k)+(-a1*(W3*C(jctrl,1)+W4*C(jctrl,2)+W3*C(jctrl,3)-a1*C(jctrl-1,1)-a2*C(jctrl-1,2)-a1*C(jctrl-1,3)+k*k*q(X(1),jctrl*k))/W1);
    for ictrl=2:elementNum+2
        B(ictrl)=W3*C(jctrl,ictrl-1)+W4*C(jctrl,ictrl)+W3*C(jctrl,ictrl+1)-a1*C(jctrl-1,ictrl-1)-a2*C(jctrl-1,ictrl)-a1*C(jctrl-1,ictrl+1)+k*k*q(X(ictrl-1),jctrl*k);
    end
    B(elementNum+3)=f2(jctrl*k)+(-a1*(W3*C(jctrl,elementNum+1)+W4*C(jctrl,elementNum+2)+W3*C(jctrl,elementNum+3)-a1*C(jctrl-1,elementNum+1)-a2*C(jctrl-1,elementNum+2)-a1*C(jctrl-1,elementNum+3)+k*k*q(X(elementNum+1),jctrl*k))/W1);
    C(jctrl+1,:)=(A\B)';
end
%% 计算U矩阵
U=zeros(timeNum+1,elementNum+1);
for jctrl=1:timeNum+1
    for ictrl=1:elementNum+1
        U(jctrl,ictrl)=a1*C(jctrl,ictrl)+a2*C(jctrl,ictrl+1)+a1*C(jctrl,ictrl+2);
    end
end
%% 后处理
plot(X,U(timeNum+1,:),X,U_EXC)
