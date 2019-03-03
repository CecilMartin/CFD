%Fromm test
%Verify that Fromm's is correctly implemented
clear
close all
clc
a=1;
L=1;
limiter=0;
rho_fun=@(x)(exp(-((x-L/2).^2)/(L/2/8)^2));
CFL=.1;

nt_num=5;
nt_base=7;
nt_all=2.^(linspace(nt_base,nt_base+nt_num-1,nt_num));

% nt=L/a/dt;
err=zeros(1,nt_num-1);
for i=1:nt_num
    nt=nt_all(i);
    dt=1/nt;
   
    dx=dt*a/CFL;
    nx=floor(L/dx);
    dx=L/nx;
    x=linspace(0,L,nx+1)';
    x=x(1:end-1);
    x_mid=x+dx/2;
    rho_fv=rho_fun(x_mid);
    for j=1:nt
        rho_fv=rho_fv+dt*Fromm(rho_fv, dt, dx, ones(nx,1)*a,limiter);
    end
    rho_real=rho_fun(x_mid);
    err(i)=norm(rho_real-rho_fv)/norm(rho_real);
end
loglog(1./nt_all,err,'-o')
hold on
% plot(nt_all, err())