function y=OneOverEhkl(Sij,psi,phi)
% Sij is compliance matrix 6x6
% psi is angle between the measured direction and [001]
% phi is angle between [001] and the projection of the measured direction to (001)

S=GetFullComplianceTensor(Sij);

y=zeros(size(psi));
for n=1:numel(psi)
    m=[sind(psi(n))*cosd(phi(n)); sind(psi(n))*sind(phi(n)); cosd(psi(n))];
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    y(n)=y(n)+S(i,j,k,l)*m(i)*m(j)*m(k)*m(l);
                end
            end
        end
    end
end