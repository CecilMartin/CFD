%Zhe Chen
%Final Project for CFD
%main script
clear;
close all;
clc;
addpath('./CodeReport/code/')

h=2;
a=0.656;
mg=1.27e-15*9.8;
dx=h/8;
M=1;
L=1024;
L_ic=12*h;%32sigma
kT=(h-a)*1e-6*mg;
T=kT/1.38064852e-23
omega=16;
S=6*a^3*omega;
Amp=M/L_ic;
%!!!!!!!!!
R=L/2;
limiter=1;
eps_f=0.025*Amp;
eps_a=.01;
V0=Amp*pi*S/3;

K=@(x)4*S*h/3*x.^2./((x.^2+4*h^2).^2);
K_tilde=@(x) K(x).*(abs(x)<R);
% rho_fun=@(x)(x<L_ic).*(x>=0).*(Amp*exp(-((x-L_ic/2).^2)/(L_ic/2/8)^2));
rho_fun=@(x)(x<L_ic).*(x>=0)*Amp;

fprintf('Truncation length of K is R/L=%.3e.\n',R/L);
if R/L>0.5
    warning('R should be less than 1/2L');
end

nx=L/dx;

fprintf('dx/h=%.3e.\n',dx/h);

x=linspace(0,L,nx+1)';
x=x(1:end-1);
% x_padd=[x-L;x];
x_mid=x+dx/2;

% rho_fv=(x_mid<L_ic)*(0.1);
% rho_fv=(x_mid<L_ic).*(0.1*exp(-((x_mid-L_ic/2).^2)/(L_ic/2)^2));
% rho_fv=(x_mid<L_ic).*(0.1*(sin(x*2*pi/L_ic-pi/2)+1));
rho_fv=rho_fun(x_mid);
rho_grid=(rho_fv+circshift(rho_fv,1))/2;

K_grid=zeros(nx,1);
K_grid(1:nx/2)=K_tilde(x(1:nx/2));
K_grid(nx/2+1:nx)=K_tilde(x(nx/2+1:nx)-L);
% figure(1)
% plot(x,K_grid,'-.o')
K_hat=fft(K_grid);

t=0;
CFL=.4;
nt=2000;
TSCREEN=nt/4;
TSCREEN_w=nt/50;
k=0;
figure;
w=zeros(1,nt/TSCREEN_w);
t_rec=zeros(1,nt/TSCREEN_w);
for t_n=1:nt
    v=real(ifft(fft(rho_grid).*K_hat))*dx;
    dt=CFL*dx/max(v);
    t=t+dt;
    
%     rho_star=rho_fv+dt*Fromm(rho_fv, dt, dx, v ,limiter);
%     rho_half=(rho_fv+rho_star)/2;
    rho_half=rho_fv+dt/2*Fromm(rho_fv, dt/2, dx, v ,limiter);
    rho_half_grid=(circshift(rho_half,1)+rho_half)/2;
    v_half=real(ifft(fft(rho_half_grid).*K_hat))*dx;
%     v_half=cconv(rho_half_grid,K_tilde(x(x<R)),nx);
    rho_fv=rho_fv+dt*Fromm(rho_fv, dt, dx, v_half,limiter );
%     rho_fv=max(rho_fv,0);
    rho_grid=(rho_fv+circshift(rho_fv,1))/2;
    if (mod(k,TSCREEN_w)==0)
        w(k/TSCREEN_w+1)=width_of_bump_3(rho_fv,dx,eps_f,eps_a);
        t_rec(k/TSCREEN_w+1)=t;
    end
    if (mod(k,TSCREEN)==0)
        plot(x_mid/h,rho_fv/Amp);
%         hold on
%         plot(x,v_half)
%         set(gca,'XLim',[0,100*h]);
        title(num2str(t));
        drawnow
        pause(.5)
    end
    k=k+1;
end
% set(gca,'XLim',[0,100*h]);
figure
plot(t_rec*V0/L_ic,w/h,'o-')
set(gca,'YLim',[4,12])
    