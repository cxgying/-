clc;clear all;close all
%% 方程定义
beta=@(x,t)(1+(x-0.5)*(x-0.5));%定义b(x)
f=@(x,t)(x*(1-x))^(1/2);%定义源项f(x)
emppower2=-1;%定义扩散系数
a=0;b=1;L=b-a;%定义求解区域长度
%% 初边值条件之定义
%初值条件
ic=@(x)(0);
dic=@(x)(0);%初始条件的空间导数
%狄利克雷边界
dbc1=@(t)(0);
dbc2=@(t)(0);
%% 配点
elementNum=100;
X=0:L/elementNum:L;
%% 求解器参数之定义
time=2;%定义求解时间
tsp=0.01;%定义时间步长
theta=0.1;
timeNum=time/tsp;%计算时间步数
%% Solve
%% 基本系数之计算
k=tsp;h=(b-a)/elementNum;
a1=sin(h/2)*sin(h/2)/(sin(h)*sin(1.5*h));
a2=2/(1+2*cos(h));
a3=-3/(4*sin(1.5*h));
a4=3/(4*sin(1.5*h));
a5=3*(1+3*cos(h))/(16*sin(0.5*h)*sin(0.5*h)*(2*cos(0.5*h)+cos(1.5*h)));
a6=-3*cos(0.5*h)*cos(0.5*h)/(2*sin(0.5*h)*sin(0.5*h)*(1+2*cos(h)));
% 待求解变量，每一行为一个时间步的数据
C=zeros(timeNum+1,elementNum+3);
%% 初始化
A=zeros(elementNum+3);
B=zeros(elementNum+3,1);
A(1,1)=a3;A(1,3)=a4;
B(1)=dic(X(1));
for ictrl=2:elementNum+2
    A(ictrl,ictrl-1)=a1;
    A(ictrl,ictrl)=a2;
    A(ictrl,ictrl+1)=a1;
    B(ictrl)=ic(X(ictrl-1));
end
A(elementNum+3,elementNum+1)=a3;A(elementNum+3,elementNum+3)=a4;
B(elementNum+3)=dic(X(elementNum+1));
C(1,:)=(A\B)';
%% 时间步进
for jctrl=1:timeNum
    A=zeros(elementNum+3);
    A(1,1)=a1;A(1,2)=a2;A(1,3)=a1;
    B(1)=0;
    for ictrl=2:elementNum+2
        bij=beta(X(ictrl-1),jctrl*k);
        bijr=beta(X(ictrl-1),(jctrl-1)*k);
        fij=f(X(ictrl-1),jctrl*k);
        fijr=f(X(ictrl-1),(jctrl-1)*k);
        A(ictrl,ictrl-1)=bij*a1+(1-theta)*k*emppower2*a5+(1-theta)*k*bij*a1;
        A(ictrl,ictrl)=bij*a2+(1-theta)*k*emppower2*a6+(1-theta)*k*bij*a2;
        A(ictrl,ictrl+1)=bij*a1+(1-theta)*k*emppower2*a5+(1-theta)*k*bij*a1;
        B(ictrl)=C(jctrl,ictrl-1)*(bijr*a1-theta*k*emppower2*a5-theta*k*bijr*a1);
        B(ictrl)=C(jctrl,ictrl)*(bijr*a2-theta*k*emppower2*a6-theta*k*bijr*a2)+B(ictrl);
        B(ictrl)=C(jctrl,ictrl+1)*(bijr*a1-theta*k*emppower2*a5-theta*k*bijr*a1)+B(ictrl);
        B(ictrl)=(1-theta)*k*fij+theta*k*fijr+B(ictrl);
    end
    A(elementNum+3,elementNum+1)=a1;A(elementNum+3,elementNum+2)=a2;A(elementNum+3,elementNum+3)=a1;
    B(elementNum+3)=0;
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
plot(X,U(1,:),X,U(0.05/k,:),X,U(0.1/k,:),X,U(0.2/k,:),X,U(0.5/k,:),X,U(1/k,:),X,U(2/k,:))


