clc;clear all;close all;
mu=0.5;
%% 网格和求解器
L=1;n=10;deltax=L/n;x=0:L/n:L;
T=1;deltat=deltax^2*mu;m=round(T/deltat);
%% 数据容器
y=zeros(m+1,n+1);
%% 定解条件
y0=sin(pi*x);
y_0=0;y_n=0;
y(1,:)=y0;
%% solve
for j=2:(T/deltat+1)
  y(j,1)=y_0;y(j,n+1)=y_n;
  for i=2:n
    y(j,i)=y(j-1,i)+mu*(y(j-1,i+1)-2*y(j-1,i)+y(j-1,i-1));
  end
end
plot(x',y(1:round((m+1)/20):(m+1),:)')
