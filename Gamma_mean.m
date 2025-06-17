N=20;
x=1:N;
y=zeros(size(x));
for n=1:N
    h=-N:1:N;
    k=-N:1:N;
    l=-N:1:N;
    [H,K,L]=meshgrid(h,k,l);
    ind=H==0 & K==0 & L==0;
    G=(H.^2.*K.^2+K.^2.*L.^2+H.^2.*L.^2)./(H.^2+K.^2+L.^2).^2;
    G1=G(ind==0);
    y(n)=mean(G1,'all');
end

plot(x,y)