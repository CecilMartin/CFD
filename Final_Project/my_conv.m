function y=my_conv(u,v)
%n>>m
%discrete convolution
n=length(u);
m=length(v);
u_padd=[zeros(m-1,1);u;zeros(m-1,1)];
y=zeros(n,1);
for i=1:n
    y(i)=sum(u_padd(i+m-1:i+m-1+m-1).*v(1:m));
    y(i)=y(i)+sum(u_padd(i-1+m-1:-1:i).*v(2:m));
end
end