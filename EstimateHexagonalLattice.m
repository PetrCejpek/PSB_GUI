function [a,c]=EstimateHexagonalLattice(hkl,d_hkl)
a0=logspace(log10(1),log10(30),20);
c0=logspace(log10(1),log10(30),20);

ind=0;
mylist=zeros(numel(a0)*numel(c0),3);
for m=1:numel(a0)
    for n=1:numel(a0)
        ind=ind+1;
        mylist(ind,:)=[a0(m) c0(n) sum(dhkl_kernel(hkl,d_hkl,[a0(m) c0(n)]).^2)];
    end
end
[~,indmin]=min(mylist(:,3));
x0=mylist(indmin,1:2);
lpT=fsolve(@(x) dhkl_kernel(hkl,d_hkl,x),x0);
a=lpT(1);
c=lpT(2);

function ss=dhkl_kernel(hkl,d_hkl,lp)
B=Reci(lp(1),lp(1),lp(2),90,90,120);
ss=zeros(size(d_hkl));
for n=1:numel(d_hkl)
    d_th=2*pi./norm(B*hkl(n,:)');
    ss(n)=(d_hkl(n)-d_th).^2;
end