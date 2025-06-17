function C=FoldComplianceTensor_4_to_2(T)
ind=[1 1;
    2 2;
    3 3;
    2 3;
    1 3;
    1 2;];

C=zeros(6,6);
for m=1:6
    for n=1:6
        C(m,n)=T(ind(m,1),ind(n,2),ind(m,1),ind(n,2));
    end
end