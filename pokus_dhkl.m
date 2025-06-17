a=4.75;
c=12.99;
B=Reci(a,a,c,90,90,120);

Hmax=10;
Kmax=10;
Lmax=10;

% [H,K,L]=meshgrid(h,k,l);
dhkl=zeros((2*Hmax+1)*(2*Kmax+1)*(2*Lmax+1),1);
n=0;
for h=-Hmax:1:Hmax
    for k=-Kmax:1:Kmax
        for l=-Lmax:1:Lmax
            n=n+1;
            dhkl(n)=norm(B*[h; k; l])/(2*pi);
        end
    end
end

plot(dhkl)