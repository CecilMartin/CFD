%Fromm
function w_t = Fromm(w, dt, dx, a ,limiter)
% a=real(ifft(fft(w).*K_hat));
% a=(circshift(a,-1)+a)/2;
if limiter
    phi=@(x)max(zeros(length(x),1),min(min((1+x)/2,2*ones(length(x),1)),2*x));
    theta=(w-circshift(w,1))./(circshift(w,-1)-w);
    fluxes=a.*circshift(w,1)+...
        0.5*a.*(1-a*dt/dx).*phi(circshift(theta,1)).*(w-circshift(w,1));
    
else
    sigma=(circshift(w,-1)-circshift(w,1))/2/dx;
    fluxes = a.*(circshift(w,1)+0.5*circshift(sigma,1).*(dx-a.*dt));
end
w_t = (fluxes -circshift(fluxes,-1))/dx;
end
