%Zhe Chen
%Final Project for CFD
%Verify convergence order of our spatial-temporal scheme
clear;
close all;
clc;
h=1;
a=0.656;
mg=1.27e-15*9.8;
CFL=.1;
L=128;
L_ic=8;%32sigma
kT=(h-a)*1e-6*mg;
T=kT/1.38064852e-23
omega=100;
S=6*a^3*omega;
Amp=1;
%!!!!!!!!!
R=L*2/3;
limiter=1;

K=@(x)4*S*h/3*x.^2./((x.^2+4*h^2).^2);
K_tilde=@(x) K(x).*(abs(x)<R);
rho_fun=@(x)(x<L_ic).*(x>=0).*(Amp*exp(-((x-L_ic/2).^2)/(L_ic/2/8)^2));
% rho_fun=@(x)(x<L_ic).*(x>=0)*Amp
% rho_fun=@(x)(x<L_ic).*(x>=0).*(Amp*(sin(x*2*pi/L_ic-pi/2)+1));

fprintf('Truncation length of K is R/L=%.3e.\n',R/L);
if R/L>0.5
    warning('R should be less than 1/2L');
end


T=.01;
nt_base=5;
nt_num=6;
nt_all=2.^linspace(nt_base,nt_base+nt_num-1,nt_num);

err=zeros(1,nt_num-1);
for j=nt_num:-1:1
    
    nt=nt_all(j);
    dt=T/nt;
    
    dx=h/4*2^(-j+1);
    nx=floor(L/dx);
    dx=L/nx;
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
    CFL=.1;
    TSCREEN=1;
    k=1;
    
    for t_n=1:nt
        v=real(ifft(fft(rho_grid).*K_hat))*dx;
%         dt=CFL*dx/max(v);
        t=t+dt;
        if dt*max(v)/dx>.5
        end
%         rho_star=rho_fv+dt*Fromm(rho_fv, dt, dx, v ,limiter);
%         rho_half=(rho_fv+rho_star)/2;
%         rho_half_grid=(circshift(rho_half,1)+rho_half)/2;
        rho_half=rho_fv+dt/2*Fromm(rho_fv, dt, dx, v/2 ,limiter);
        rho_half_grid=(circshift(rho_half,1)+rho_half)/2;
        
        v_half=real(ifft(fft(rho_half_grid).*K_hat))*dx;
        %     v_half=cconv(rho_half_grid,K_tilde(x(x<R)),nx);
%         v_half=ones(nx,1);
        rho_fv=rho_fv+dt*Fromm(rho_fv, dt, dx, v_half,limiter );
%         rho_fv=rho_star;
        %     rho_fv=max(rho_fv,0);
        rho_grid=(rho_fv+circshift(rho_fv,1))/2;
        
%         if (mod(k,TSCREEN)==0)
%             plot(x_mid,rho_fv);
% %             plot(x,v_half)
%             %             set(gca,'XLim',[0,10*h]);
%             title(num2str(t));
%             drawnow
% %             pause(1)
%         end
        k=k+1;
        
    end
    t
    
    if j~=nt_num
        err(j)=norm(rho_fv-rho_finer)/norm(rho_finer);
    end
%     rho_finer=rho_fv;
    rho_finer=Coarsen(rho_fv,2);
end
figure;
log(err(end)/err(end-1))/log(nt_all(end-1)/nt_all(end-2))
loglog(nt_all(1:end-1),err,'o','DisplayName','Relative error')
hold on
plot(nt_all(1:end-1),err(end)*(nt_all(1:end-1)/nt_all(end-1)).^-2,'DisplayName','2nd order')
legend show