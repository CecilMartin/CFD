function w=width_of_bump_3(rho,dx,eps_f,eps_a)
%Compute width of front shock
n=length(rho);
a=rho<eps_f;
for i=n:-1:2
    a(1:i-1)=a(1:i-1)*a(i);
end
fore=find(a,1);
% rho_x=(circshift(rho,-1)-circshift(rho,1))/2/dx;
flag=0;
for i=n:-1:1
    if flag==0 
        if (rho(i))<eps_f
            continue
        else
            flag=1;
        end
    elseif flag==1
        if rho(i)>eps_f
            continue
        else
            flag=2;
            aft=i;
            break;
        end
    end
end
if flag==1
    aft=1;
end
w=(fore-aft)*dx;
end