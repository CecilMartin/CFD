%Regularization Test
%Zhe Chen
clear;
close all;
clc;
addpath('./CodeReport/code/')
%------------------------------------------------
%parameters
h=2;
a=0.656;
mg=1.27e-15*9.8;
% dx=h/8;
dx_all=h.*(2.^(-5:-1));
M=1;
L=512;
L_ic=12*h;%32sigma
kT=(h-a)*1e-6*mg;
T=kT/1.38064852e-23
omega=16;
S=6*a^3*omega
Amp=M/L_ic;
%!!!!!!!!!
R=L/2;
limiter=1;
eps_f=0.025*Amp;
eps_a=.01;
V0=Amp*pi*S/3;
BDS_flag=1;
sigma=h*2;
% 
% sigma_all=2^(-2:)
% K=@(x)4*S*h/3*x.^2./((x.^2+4*h^2).^2);
K=@(x)(4*(4*h^4+3*h^2*x.^2)./((4*h^2+x.^2).^2)+log(1+4*h^2./(x.^2)))/8/pi;
K_tilde=@(x) K(x).*(abs(x)<R);
% rho_fun=@(x)(x<L_ic).*(x>=0).*(Amp*exp(-((x-L_ic/2).^2)/(L_ic/2/8)^2));
rho_fun=@(x)(x<L_ic).*(x>=0)*Amp;
gauss_kernal=@(x)(1/(2*pi*sigma^2))^(1/2)*exp(-x.^2/(2*sigma^2));
err=zeros(length(dx_all)-1,1);
fprintf('Truncation length of K is R/L=%.3e.\n',R/L);
if R/L>0.5
    warning('R should be less than 1/2L');
end
%-----------------------------------------------
%Start
figure
hold on
for i=length(dx_all):-1:1
    dx=dx_all(i);
    nx=L/dx;
    fprintf('dx/h=%.3e.\n',dx/h);
    x=linspace(0,L,nx+1)';
    x=x(1:end-1);
    x_mid=x+dx/2;
    K_mid=zeros(nx,1);
    K_mid(1:nx/2)=K_tilde(x_mid(1:nx/2));
    K_mid(nx/2+1:nx)=K_tilde(x_mid(nx/2+1:nx)-L);
    gauss=zeros(nx,1);
    gauss(1:nx/2)=gauss_kernal(x_mid(1:nx/2));
    gauss(nx/2+1:nx)=gauss_kernal(x_mid(nx/2+1:nx)-L);
    K_hat=fft(K_mid).*fft(gauss);
    K_sigma=real(ifft(K_hat))*dx;
    plot(x_mid,K_sigma,'DisplayName',['dx/h=1/',num2str(h/dx)]);
    if i~=length(dx_all)
        err(i)=norm(Coarsen(K_sigma,2)-K_sigma_old,2)/norm(Coarsen(K_sigma,2),2);
    end
    K_sigma_old=K_sigma;
end
legend show
set(gca,'FontSize',16);
figure
hold on
% plot(x_mid,K_mid)
% plot(x_mid,K_sigma,'o--');
plot(dx_all(2:end),err,'o-','DisplayName','Relative Error')
plot(dx_all(2:end),dx_all(2:end).^2*err(1)/dx_all(2)^2,'--','DisplayName','2nd order');
plot(dx_all(2:end),dx_all(2:end).^1*err(1)/dx_all(2)^1,'-.','DisplayName','1st order');
set(gca,'YScale','log','XScale','log')
set(gca,'fontsize',16)
xlabel('dx');ylabel('Relative norm 2 Error')
title('Emphirical Relative Error')
legend show
