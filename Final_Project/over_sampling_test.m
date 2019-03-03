%V test
%Zhe Chen
%Final Project for CFD
clear;
close all;
clc;

h=1;
a=0.656;
mg=1.27e-15*9.8;

dx=h/2;
L_ic=16;
R=32; %R should be less than 1/2L
kT=(h-a)*1e-6*mg;
T=kT/1.38064852e-23
omega=100;
S=6*a^3*omega;


Amp=1;
coef_a=Amp*4*h*S/3;
coef_b=2*h;
L_all=2.^([6,7,8,9,10]);
err=zeros(1,length(L_all));


K=@(x)4*S*h/3*x.^2./((x.^2+4*h^2).^2);
K_tilde=@(x) K(x).*(abs(x)<R);
rho_fun=@(x)(x<L_ic).*(x>=0).*(Amp*exp(-((x-L_ic/2).^2)/(L_ic/2/8)^2));
v_true_fun=@(x)integral(@(y)rho_fun(x-y).*K(y),-R,R,'AbsTol',1e-12,'RelTol',1e-12);
% v_true_fun=@(x)coef_a/2*((atan(x/coef_b)-atan((x-L_ic)/coef_b))/coef_b-...
%     x./(coef_b^2+x.^2)+(x-L_ic)./(coef_b^2+(x-L_ic).^2));

% rho_fun=@(x)(x<L_ic).*(x>=0).*(Amp*(sin(x*2*pi/L_ic-pi/2)+1));

% rho_fun_1=@(x)(x<L_ic).*(x>=0).*(Amp*exp(-((x-L_ic/2).^2)/(L_ic/2)^2));
% rho_fun=@(x)rho_fun_1(x)-rho_fun_1(0);
    

for i=1:length(L_all)
    L=L_all(i);
    fprintf('Truncation length of K is R/L=%.3e.\n',R/L);
    if R/L>0.5
        warning('R should be less than 1/2L');
    end
    nx=L/dx;
    fprintf('dx/h=%.3e.\n',dx/h);
    
    x=linspace(0,L,nx+1)';
    x=x(1:end-1);
    % x_padd=[x-L;x];
    % x_mid=x+h/2;
    %     rho=(x<L_ic)*Amp;
    % rho_fv=(x_mid<L_ic).*(0.1*exp(-((x_mid-L_ic/2).^2)/(L_ic/2)^2));
    
    rho=rho_fun(x);
    
    % figure(1)
    % plot(x,K_tilde(x),'-.o')
    K_grid=zeros(nx,1);
    K_grid(1:nx/2)=K_tilde(x(1:nx/2));
    K_grid(nx/2+1:nx)=K_tilde(x(nx/2+1:nx)-L);
    K_hat=fft(K_grid);
    
    
    
    v=real(ifft(fft(rho).*K_hat))*dx;
    % v(nx/2+1:nx)=0;
    % v=conv(rho,K_tilde(x(x<R)),'same')*dx;
    % v=my_conv(rho,K_tilde(x(x<R)))*dx;
    
    
    % figure
    % hold on
    % set(gca,'XLim',[-100*h,100*h]);
    index_eval=(x<R);
    % plot(x(index_eval),v(index_eval),'o');
    % fplot(v_true_fun,[-100*h,100*h]);
    % plot(x,v_true_fun(x));
    v_true=zeros(sum(index_eval),1);
    for j=1:sum(index_eval)
        v_true(j)=v_true_fun(x(j));
    end
    err(i)=norm(v(index_eval)-v_true)/norm(v_true);
    % plot(x(index_eval),v(index_eval)-v_true_fun(x(index_eval)))
    
end
figure;
hold on
% plot(L_all/R,err(1)*(L_all/L_all(1)).^-1,'DisplayName','order 1');
% plot(L_all/R,err(1)*(L_all/L_all(1)).^-2,'DisplayName','order 2');

plot(L_all/R,err,'o--','DisplayName','Relative Error')
set(gca,'XScale','log','YScale','log')
legend show


