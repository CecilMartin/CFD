clear all
close all
L = 1;
CFL = .2;%0.2;
tfinal = 1;


%----------------- initial condition ---------------

%--------------------------------------
K = 1;
limit = 1;
type = 'BDS';
limiter=0;
res = 2.^[5:9];%[5:8];
error = zeros(4,length(res));
u_unlimit = cell(1,length(res));
% u_limit = cell(1,length(res));
% u_unlimit2 = cell(1,length(res));
% u_limit2 = cell(1,length(res));
u_exact = cell(1,length(res));

c = 2;
fa  = @(x) sin(2*pi*x)+c; % advection function
sol =  @(x,t) cos(2*pi*(x+t)); % MMS
sol_aver = @(x,t) sin(2*pi*(x+t))/2/pi;

S = @(x,t) 2*pi*cos(2*pi*(2*x+t))-(2+2*c)*pi*sin(2*pi*(x+t));
Sint = @(x,t) -cos(2*pi*(2*x+t))/4/pi+(1+c)/2/pi*sin(2*pi*(x+t));
figure;
for j=1:length(res)
    
    N = res(j);
    dx = L/N;
    
    x = (0:N-1)*dx+1*dx/2; % centers of cells
    xf = x+dx/2; % faces of cells
    %----------------------- initialize ----------------------------
    
    af = fa(xf'); %face values
    a = fa(x');
    
    %     %---test---
    %     sol = S;
    %     af = 0*af;
    %     a = 0*a;
    %
    %     %----------
    
    
    %u = sol(x',0) + 1/24*( sol(x([N,1:N-1])',0)+sol(x([2:N,1])',0)-2*sol(x',0)) ;%
    u = ( sol_aver(xf',0)-sol_aver(xf([N,1:N-1])'  ,0))/dx; %f(x');
    
    u_unlimit{j} = u;
    %     u_limit{j} = u;
    %     u_unlimit2{j} = u;
    %     u_limit2{j} = u;
    %----------------- end of initialize ---------------------------
    
    %-------------------discretization ---------------
    
    
    
    r = CFL;
    dt = r*dx;
    %tfinal = 1*dt ; % test for one step error
    nt = round(tfinal/dt);
    % ------------------ end of advection coefficient ------------------
    
    %u_exact = sol(x',tfinal );
    %u_exact{j} = sol(x',tfinal) + 1/24*( sol(x([N,1:N-1])',tfinal)+sol(x([2:N,1])',tfinal)-2*sol(x',tfinal)) ;
    u_exact{j} = ( sol_aver(xf',tfinal)-sol_aver(xf([N,1:N-1])'  ,tfinal))/dx; %f(x');
    t0 = 0;
    %figure;
    for i = 1:nt
        
        
        St = Sint(xf',t0+dt)-Sint(xf',t0)-Sint( xf([N,1:N-1])' ,t0+dt)+Sint(xf([N,1:N-1])',t0) ;
        St = St/dt/dx;
        
        
        t0 = t0 + dt;
        %St = St*0;
        rho_half=u_unlimit{j}+dt/2*Fromm(u_unlimit{j}, dt/2, dx, af ,limiter)+St*dt/2;
        %         rho_half_grid=(circshift(rho_half,1)+rho_half)/2;
        %         v_half=real(ifft(fft(rho_half_grid).*K_hat))*dx;
        %     v_half=cconv(rho_half_grid,K_tilde(x(x<R)),nx);
        
        u_unlimit{j}=u_unlimit{j}+dt*Fromm(rho_half, dt, dx,af,limiter )+St*dt;
        
        %         u_unlimit{j} = temporal_integrator(u_unlimit{j},af,dx,dt,1,0,type,St);
        %flux_exact = fa(x'+dx/2).*(1/2/pi/dt).*( sin(2*pi*(x'+dx/2+t0+dt)) -  sin(2*pi*(x'+dx/2+t0)) );
        %error_flux(j) = norm(flux-flux_exact)/norm(flux_exact);
        %aver_exact = (1/2/pi/dt).*( sin(2*pi*(x'+dx/2+t0+dt)) -  sin(2*pi*(x'+dx/2+t0)) );
        %error_aver(j) = norm(s_aver-aver_exact)/norm(aver_exact);
        % u_unlimit{j} = temporal_integrator_test(u_unlimit{j},af,dx,dt,1,0,type,St);
        %         u_limit{j} = temporal_integrator(u_limit{j},af,dx,dt,1,1,'BDS',St);
        %         u_unlimit2{j} = temporal_integrator(u_unlimit2{j},af,dx,dt,2,0,type,St);
        %         u_limit2{j} = temporal_integrator(u_limit2{j},af,dx,dt,2,1,'M',St);
        
    end
    error(1,j) = norm(u_exact{j}-u_unlimit{j},1)/norm(u_exact{j},1);
end

%%
%stop
close all
line_spec = {'o-.','*:','+-','s--'};

% 
% figure;
% %subplot(1,3,1)
% for i = 1:1
%     loglog(res,error(i,:),line_spec{i},'linewidth',2),hold on,
% end
% hold on
% xlabel('N')
% ylabel('\epsilon')
% title('global relative error')
% set(gca,'fontsize',16)
% yl = ylim;
% 
% xl = [res(1)/2,res(end)*2];
% 
% 
% intp(1,:) = polyfit(log(1./res(2:end)),log(error(1,2:end)),1);
% intp(2,:) = polyfit(log(1./res(2:end)),log(error(2,2:end)),1);
% intp(3,:) = polyfit(log(1./res(2:end)),log(error(3,2:end)),1);
% intp(4,:) = polyfit(log(1./res(2:end)),log(error(4,2:end)),1);
% 
% hold on;
% 
% loglog(res,  (res./res(1)).^-1*error(1) ,'k--','linewidth',1)
% loglog(res,  (res./res(1)).^-2*error(1) ,'k--','linewidth',1)
% loglog(res,  (res./res(1)).^-3*error(1) ,'k--','linewidth',1)
% % loglog(xl,  exp(intp(2,2))* (1./xl).^2/xl(1) ,'k-.','linewidth',1)
% % loglog(xl,  exp(intp(2,2))* (1./xl).^3*xl(1) ,'k-','linewidth',1)
% %         legend('unlimited, linear','limited, linear','unlimited, quadratic','limited, quadratic','slope=-1','slope=-2','slope=-3')
% 
% 
% index = find(res==64);
% N = res(index);
% dx = L/N;
% dt = r*dx;
% nt = round(tfinal/dt);
% x = (0:N-1)*dx+dx/2;
% 
% a = fa(x');


figure;

for index = 2:length(res)
    N = res(index);
    dx = L/N;
    x = (0:N-1)*dx+dx/2;
    plot(x,u_unlimit{index}-u_exact{index},line_spec{index-1}),hold on,
    
end
legend('N=64','N=128','N=256','N=512')
xlim([x(1),x(end)])
%ylim([-0.2,1.2])
title(error)
yl=ylabel('$s_h-s_{exact}$','interpreter','latex')
title('error with linear reconstruction')
xlabel('x')
set(gca,'fontsize',16)

% subplot(1,2,2)
% for index = 2:length(res)
%     N = res(index);
%     dx = L/N;
%     x = (0:N-1)*dx+dx/2;
%     plot(x,u_unlimit2{index}-u_exact{index},line_spec{index-1}),hold on,
% end
legend('N=64','N=128','N=256','N=512')
xlim([x(1),x(end)])
%ylim([-0.2,1.2])
title(error)
yl=ylabel('$s_h-s_{exact}$','interpreter','latex')
title('error with quadratic reconstruction')
xlabel('x')
set(gca,'fontsize',16)



%% empirical order estimation


ratio = zeros(4,length(res)-1);
for i = 1:length(res)-1
    ratio(:,i) = log( error(:,i+1)./error(:,i))/log(1/2);
end
ratio





