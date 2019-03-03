function w=width_of_bump_2(rho,dx,eps_f,eps_a)
n=length(rho);
a=rho<eps_f;
for i=n:-1:2
    a(1:i-1)=a(1:i-1)*a(i);
end
fore=find(a,1);
rho_x=(circshift(rho,-1)-circshift(rho,1))/2/dx;
flag=0;
for i=n:-1:1
    if flag==0 
        if abs(rho(i))<eps_f
            continue
        else
            flag=1;
        end
    elseif flag==1
        if rho_x(i)<0
            continue
        else
            if abs(rho_x(i-1))<eps_a
                flag_2=1;
            else
                flag_2=0;
                flag=3;
                continue;
            end
            flag=2;
        end
    elseif flag==2
        if flag_2==1
            if abs(rho_x(i))<eps_a
                continue
            else
                flag=3;
            end
        else
            if flag~=3
                error('tmd')
            end
        end
    elseif flag==3
        if abs(rho_x(i))<eps_a
            aft=i;
            flag=4;
            break
        end
    end
end
if flag==3
    aft=1;
end
w=(fore-aft)*dx;
end