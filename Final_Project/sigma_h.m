%Zhe Chen
%Final Project for CFD
%main script
clear;
close all;
clc;
addpath('./CodeReport/code/')
%------------------------------------------------
%parameters
h_all=2.^(2:3);
a=0.656;
mg=1.27e-15*9.8;
dx=1/8;
M=1;
L=256;
L_ic=24;%32sigma
% kT=(h-a)*1e-6*mg;
% T=kT/1.38064852e-23
omega=16;
S=6*a^3*omega
Amp=M/L_ic;
%!!!!!!!!!
R=L/2;
limiter=1;
eps_f=0.025*Amp;
eps_a=3;
V0=Amp*pi*S/3;
BDS_flag=1;
sigma_all=2.^(0:2);
width_rec=zeros(length(h_all),length(sigma_all));
if R/L>0.5
    warning('R should be less than 1/2L');
end
fprintf('Truncation length of K is R/L=%.3e.\n',R/L);
nx=L/dx;
% fprintf('dx/h=%.3e.\n',dx/h);
rho_fun=@(x)(x<L_ic).*(x>=0)*Amp;
for i=1:length(h_all)
    i
    for j=1:length(sigma_all)
        j
        h=h_all(i);
        sigma=sigma_all(j);
        K=@(x)(4*(4*h^4+3*h^2*x.^2)./((4*h^2+x.^2).^2)+log(1+4*h^2./(x.^2)))/8/pi;
        K_tilde=@(x) K(x).*(abs(x)<R);
        x=linspace(0,L,nx+1)';
        x=x(1:end-1);
        x_mid=x+dx/2;
        
        % rho_fv=(x_mid<L_ic)*(0.1);
        % rho_fv=(x_mid<L_ic).*(0.1*exp(-((x_mid-L_ic/2).^2)/(L_ic/2)^2));
        % rho_fv=(x_mid<L_ic).*(0.1*(sin(x*2*pi/L_ic-pi/2)+1));
        rho_fv=rho_fun(x_mid);
        % rho_grid=(rho_fv+circshift(rho_fv,1))/2;
        gauss_kernal=@(x)(1/(2*pi*sigma^2))^(1/2)*exp(-x.^2/(2*sigma^2));
        
        K_mid=zeros(nx,1);
        K_mid(1:nx/2)=K_tilde(x_mid(1:nx/2));
        K_mid(nx/2+1:nx)=K_tilde(x_mid(nx/2+1:nx)-L);
        
        gauss=zeros(nx,1);
        gauss(1:nx/2)=gauss_kernal(x_mid(1:nx/2));
        gauss(nx/2+1:nx)=gauss_kernal(x_mid(nx/2+1:nx)-L);
        K_hat=fft(K_mid).*fft(gauss)*dx;
        % figure
        % hold on
        % plot(x_mid,real(ifft(K_hat)))
        % plot(x_mid,K_mid)
        t=0;
        CFL=.4;
        nt=2000;
        TSCREEN=nt/4;
        TSCREEN_w=nt/50;
        k=0;
        dt=0.4;
%         figure;
        w=zeros(1,nt/TSCREEN_w);
        t_rec=zeros(1,nt/TSCREEN_w);
        max_rec=zeros(1,nt/TSCREEN_w);
        for t_n=1:nt
            v=real(ifft(fft(rho_fv).*K_hat))*dx;
            v=(v+circshift(v,1))/2;
%             dt=CFL*dx/max(v);
            t=t+dt;
%             if dt>CFL*dx/max(v)
%                 error('cnm')
%             end
            if BDS_flag
                rho_fv=temporal_integrator(rho_fv,v,dx,dt,2,1,'PPM2',0);
            else
                %     rho_star=rho_fv+dt*Fromm(rho_fv, dt, dx, v ,limiter);
                %     rho_half=(rho_fv+rho_star)/2;
                rho_half=rho_fv+dt/2*Fromm(rho_fv, dt/2, dx, v ,limiter);
                rho_half_grid=(circshift(rho_half,1)+rho_half)/2;
                v_half=real(ifft(fft(rho_half_grid).*K_hat))*dx;
                %     v_half=cconv(rho_half_grid,K_tilde(x(x<R)),nx);
                rho_fv=rho_fv+dt*Fromm(rho_fv, dt, dx, v_half,limiter );
                %     rho_fv=max(rho_fv,0);
                
                %
            end
            %     rho_grid=(rho_fv+circshift(rho_fv,1))/2;
            
            if (mod(k,TSCREEN_w)==0)
                w(k/TSCREEN_w+1)=width_of_bump(rho_fv,dx,eps_f,eps_a*Amp,2);
                t_rec(k/TSCREEN_w+1)=t;
                max_rec(k/TSCREEN_w+1)=max(rho_fv);
            end
%             if (mod(k,TSCREEN)==0)
%                 plot(x_mid/h,rho_fv/Amp);
%                 %         hold on
%                 %         plot(x,v_half)
%                 %         set(gca,'XLim',[0,100*h]);
%                 title(num2str(t));
%                 drawnow
%                 pause(.5)
%             end
            k=k+1;
            end0-
        width_rec(i,j)=mean(w(end/5*4:end));
    end
end

%%
figure
hold on
linspec=['o','+']
for i=1:length(h_all)
   plot(sigma_all,width_rec(i,:),linspec(i),'DisplayName',['h=',num2str(h_all(i))])
end
legend show
set(gca,'FontSize',16)
xlabel('\sigma')
ylabel('Width of front shock')
        
