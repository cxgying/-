clc;clear all;close all
%% 参数之预定义
x_alpha_G=0.2;
r_alpha_G=0.5;
mu=20;
a_G=-0.1;
R_omega=0.3;
rho_alpha=1.225;
%% Theodosren函数句柄之定义
C=@(k)(besselh(1,2,k)/(besselh(1,2,k)+1i*besselh(0,2,k)));
%% 循环参数
k_set=0.05:0.002:10;
ictrl=1;dataA=zeros(10,2);
%% 计算广义特征值问题
for k=k_set
    % 内部参数之计算
    L_h=(1-1i*2*C(k)/k);
    M_h=0.5;
    L_alpha=0.5-1i*(1+2*C(k))/k-2*C(k)/(k*k);
    M_alpha=0.375-1i/k;
    % 矩阵之定义
    M=[1 x_alpha_G;x_alpha_G r_alpha_G^2];
    K=[R_omega^2 0;0 r_alpha_G^2];
    A=[L_h L_alpha-(0.5+a_G)*L_h;M_h-(0.5+a_G)*L_h M_alpha-(0.5+a_G)*(L_alpha+M_h)+L_h*(0.5+a_G)^2]/mu;
    % 广义特征值之计算
    A=A+M;
    A_poly=K(1,1)*K(2,2);B_poly=-K(1,1)*A(2,2)-A(1,1)*K(2,2);C_poly=A(1,1)*A(2,2)-(A(1,2)-K(1,2))*(A(2,1)-K(2,1));
    H(2)=(-B_poly-sqrt(B_poly*B_poly-4*A_poly*C_poly))/(2*A_poly);
    H(1)=(-B_poly+sqrt(B_poly*B_poly-4*A_poly*C_poly))/(2*A_poly);
    %%后处理系数之计算
    % 扭转分支
    lamda_ReA=real(H(2));lamda_ImA=imag(H(2));
    V_G=1/(k*sqrt(lamda_ReA));g=lamda_ImA/lamda_ReA;omega_G=1/sqrt(lamda_ReA);
    dataA(ictrl,:)=[V_G g];dataC(ictrl,:)=[V_G omega_G];
    lamda_ReB=real(H(1));lamda_ImB=imag(H(1));
    V_G=1/(k*sqrt(lamda_ReB));g=lamda_ImB/lamda_ReB;omega_G=1/sqrt(lamda_ReB);
    dataB(ictrl,:)=[V_G g];dataD(ictrl,:)=[V_G omega_G];
    ictrl=ictrl+1;
end
%% 寻找函数零点
f=@(zeros)(spline(dataA(:,1),dataA(:,2),zeros));
fzero(f,6)
%%
%% 后处理
figure()%画第一张图
plot(dataA(:,1),dataA(:,2),dataB(:,1),dataB(:,2))
figure()%画第二张图
plot(dataC(:,1),dataC(:,2),dataD(:,1),dataD(:,2))
