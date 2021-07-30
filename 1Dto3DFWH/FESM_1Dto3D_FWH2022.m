%% Author: Wang Lican; Lab: XTERs
clear all; close all;clc
M1=0.0;M2=0.0;M3=0.0;  
%% 单极子声源
x1=[0];x2=[0];x3=[0];                    	% 声源位置
x_o=linspace(10,50,10);y_o=0*x_o;z_o=0*x_o;                        % 观测点坐标
%x_o=100;y_o=0*x_o;z_o=0*x_o;                        % 观测点坐标
rho_0=1.225;c0=340;omega=10*pi;Amp=1;       % 基本参数
dt=2*pi/omega/50;   t=0:dt:10;              % t=(1:5000)*dt;             % 时间
wlc_fmin=50;wlc_fmax=54;
dim=3;
if dim==1
    %% 1D 点和等效源面
    [x_scat, y_scat, z_scat, dS_scat, nx_scat, ny_scat, nz_scat] = datasurface1D(0,0,0,1);    
    [x_esm, y_esm, z_esm, dS_esm, nx_esm, ny_esm, nz_esm] = datasurface1D(0,0,0,0);    
elseif dim==2
    %% 2D 圆和等效源面
    [x_scat, y_scat, z_scat, dS_scat, nx_scat, ny_scat, nz_scat] = datasurface2D(0,0,0,0.2,500);    
    [x_esm, y_esm, z_esm, dS_esm, nx_esm, ny_esm, nz_esm] = datasurface2D(0,0,0,0.8*0.2,160);
elseif dim==3
    %% 3D 球面和等效源面
    [x_scat, y_scat, z_scat, dS_scat, nx_scat, ny_scat, nz_scat] = datasurface(0,0,0,1    ,2);    
    [x_esm , y_esm , z_esm , dS_esm , nx_esm , ny_esm , nz_esm ] = datasurface(0,0,0,1*0.8,1);    
end
%% 频域参数
fs_data=1/dt;   % 采样频率
N=length(t);NH=length(x_scat);Nj=length(x_esm);Nt=length(t); % 信号长度
f_data=(0:N-1)*fs_data/N; % 离散频率
%% 壁面边界信号 (用FW-H方程得到壁面上声源点信息)
P_ESM=zeros(length(x_scat),length(t));              % 存储三个方向的压力梯度
gamma=1/sqrt(1-M1^2-M2^2-M3^2);                     % M1 M2 M3表示文献[1]中的M∞ 
for j1=1:length(x_scat)
    up1=zeros(1,length(t));up2=zeros(1,length(t));up3=zeros(1,length(t));pGw_pt=zeros(1,length(t));pp=zeros(1,length(t));
    x=x_scat(j1)-x1;y=y_scat(j1)-x2;z=z_scat(j1)-x3;
    r=sqrt(x.^2+y.^2+z.^2);                            % 声源点和积分面上点的直线距离
    R_s=1/gamma*sqrt(r^2+gamma^2*(M1*x+M2*y+M3*z)^2);  % R_s表示R*（s表示star的意思）
    R=gamma^2*(R_s-M1*x-M2*y-M3*z);                    
    pr_px=x/r;pr_py=y/r;pr_pz=z/r;

    pRs_px=1/gamma^2/R_s*(r*pr_px+gamma^2*(M1*x+M2*y+M3*z)*M1);    % pRs_px 表示R*对x的偏导数（p为partial）
    pRs_py=1/gamma^2/R_s*(r*pr_py+gamma^2*(M1*x+M2*y+M3*z)*M2);    % pRs_py 表示R*对y的偏导数（p为partial）
    pRs_pz=1/gamma^2/R_s*(r*pr_pz+gamma^2*(M1*x+M2*y+M3*z)*M3);    % pRs_pz 表示R*对y的偏导数（p为partial）
    Mr=M1*x+M2*y+M3*z;
    G_w= -Amp*1i*gamma/4*(omega/c0/2/pi/R_s).^(dim/2-1)*besselh(dim/2-1,2,omega/c0*gamma^2*R_s)*exp(1i*omega/c0*gamma^2*Mr+1i*omega*t);    
    pGw_pt=pGw_pt+ 1i*omega*G_w;
    up1=up1+ ((1-dim/2)*pRs_px./R_s+1i*omega/c0*gamma^2*M1+omega/c0*gamma^2*pRs_px*(besselh(dim/2-2,2,omega/c0*gamma^2*R_s)-besselh(dim/2,2,omega/c0*gamma^2*R_s))/2/besselh(dim/2-1,2,omega/c0*gamma^2*R_s))*G_w;
    up2=up2+ ((1-dim/2)*pRs_py./R_s+1i*omega/c0*gamma^2*M2+omega/c0*gamma^2*pRs_py*(besselh(dim/2-2,2,omega/c0*gamma^2*R_s)-besselh(dim/2,2,omega/c0*gamma^2*R_s))/2/besselh(dim/2-1,2,omega/c0*gamma^2*R_s))*G_w;
    up3=up3+ ((1-dim/2)*pRs_pz./R_s+1i*omega/c0*gamma^2*M3+omega/c0*gamma^2*pRs_pz*(besselh(dim/2-2,2,omega/c0*gamma^2*R_s)-besselh(dim/2,2,omega/c0*gamma^2*R_s))/2/besselh(dim/2-1,2,omega/c0*gamma^2*R_s))*G_w;
    pp =pp-rho_0*(pGw_pt+M1*c0*up1+M2*c0*up2+M3*c0*up3);

    [pGw_px,pGw_py,pGw_pz,pGw_pxpx,pGw_pypx,pGw_pzpx,pGw_pxpy,pGw_pypy,pGw_pzpy,pGw_pxpz,pGw_pypz,pGw_pzpz] = GreenFunctionOrder2(x,y,z,M1,M2,M3,dim,omega,c0,Amp,t);
    dPdx(j1,:) = -rho_0*(1i*omega*pGw_px+M1*c0*pGw_pxpx+M2*c0*pGw_pypx+M3*c0*pGw_pzpx);
    dPdy(j1,:) = -rho_0*(1i*omega*pGw_py+M1*c0*pGw_pxpy+M2*c0*pGw_pypy+M3*c0*pGw_pzpy);
    dPdz(j1,:) = -rho_0*(1i*omega*pGw_pz+M1*c0*pGw_pxpz+M2*c0*pGw_pypz+M3*c0*pGw_pzpz);

    % 将压力梯度转为频域信号   
    up1_ESM(j1,:)=up1;
    up2_ESM(j1,:)=up2;
    up3_ESM(j1,:)=up3;    
    pp_ESM(j1,:)=pp;
    Rho_ESM(j1,:)=pp/c0^2+rho_0;
end
%% 等效源法+FW-H方程
tic
for j1=1:length(x_scat)
    P_ESM(j1,:)=fft(pp_ESM(j1,:),N);  % 不加负号，因为压力大小相等方向相同
end
pp_sum1=FWH2ESM_module(x_scat, y_scat, z_scat, x_esm, y_esm, z_esm, ...
     M1, M2, M3, dim, c0, t, P_ESM, wlc_fmin, wlc_fmax, x_o, y_o, z_o);
toc 
%% FW-H方程
tic
pp_sum2=FWH_module(x_scat,y_scat,z_scat,dS_scat,x_o,y_o,z_o,...
     nx_scat, ny_scat, nz_scat, up1_ESM,up2_ESM,up3_ESM,pp_ESM,Rho_ESM,wlc_fmin,wlc_fmax,M1,M2,M3,c0,rho_0,dt,dim);
toc
%% 等效源法+声散射法
tic
for j1=1:length(x_scat)
    dPdx_ESM(j1,:)=fft( dPdx(j1,:),N);  
    dPdy_ESM(j1,:)=fft( dPdy(j1,:),N);  
    dPdz_ESM(j1,:)=fft( dPdz(j1,:),N);      
end
[pdpdx_sum1,pdpdy_sum1,pdpdz_sum1]=SCAT2ESM_module(x_scat, y_scat, z_scat, nx_scat, ny_scat, nz_scat, x_esm, y_esm, z_esm, ...
     M1, M2, M3, dim, c0, t, dPdx_ESM,dPdy_ESM,dPdz_ESM, wlc_fmin, wlc_fmax, x_o, y_o, z_o);
toc 
%% 声散射法
tic
[pdpdx_sum2,pdpdy_sum2,pdpdz_sum2]=SCAT_module(x_scat,y_scat,z_scat,dS_scat,x_o,y_o,z_o,...
     nx_scat, ny_scat, nz_scat, up1_ESM,up2_ESM,up3_ESM,pp_ESM,Rho_ESM,wlc_fmin,wlc_fmax,M1,M2,M3,c0,rho_0,dt,dim);
toc 

source_num=10;
plot(t/(2*pi/omega),real(pp_sum1(source_num,:)),'r','linewidth',2)
hold on 
scatter(t/(2*pi/omega),real(pp_sum2(source_num,:)),'b','linewidth',2)
%% 入射波解析解（没问题）
pp = GreenFunctionOrder0(x_o(source_num)-x1,y_o(source_num)-x2,z_o(source_num)-x3,M1,M2,M3,dim,omega,c0,Amp,t,rho_0);
plot(t/(2*pi/omega),real(pp),'b--','linewidth',2)
% err=sqrt(sum((real(pp_sum1)-real(pp)).^2)/sum(real(pp).^2))*100;
xlim([20 21])
figure(2)
plot(t/(2*pi/omega),real(pdpdx_sum1(source_num,:)),'r','linewidth',2)
hold on 
scatter(t/(2*pi/omega),real(pdpdx_sum2(source_num,:)),'b','linewidth',2)
[pGw_px,pGw_py,pGw_pz,pGw_pxpx,pGw_pypx,pGw_pzpx,pGw_pxpy,pGw_pypy,pGw_pzpy,pGw_pxpz,pGw_pypz,pGw_pzpz] = GreenFunctionOrder2(x_o(source_num)-x1,y_o(source_num)-x2,z_o(source_num)-x3,M1,M2,M3,dim,omega,c0,Amp,t);
dPdx_ESM = -rho_0*(1i*omega*pGw_px+M1*c0*pGw_pxpx+M2*c0*pGw_pypx+M3*c0*pGw_pzpx);
dPdy_ESM = -rho_0*(1i*omega*pGw_py+M1*c0*pGw_pxpy+M2*c0*pGw_pypy+M3*c0*pGw_pzpy);
dPdz_ESM = -rho_0*(1i*omega*pGw_pz+M1*c0*pGw_pxpz+M2*c0*pGw_pypz+M3*c0*pGw_pzpz);
plot(t/(2*pi/omega),real(dPdx_ESM),'b--','linewidth',2)

xlim([20 21])

%% 保存数据
% var_for_plot(1,:)=t/(2*pi/omega);
% var_for_plot(2,:)=real(pp_sum1);
% var_for_plot(3,:)=t/(2*pi/omega);
% var_for_plot(4,:)=real(pp);