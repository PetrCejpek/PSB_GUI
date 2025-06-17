function odf=ODFevaluate(phi1,psi,phi2,ODF,sym)

switch ODF.TextureType
    case 'fibre'
        % distance from the line connecting or1 and or2
        or1=ODF.parameters.or1;
        or2=ODF.parameters.or2;
        sigma=ODF.parameters.sigma;
        d=zeros(size(phi1));
        for n=1:numel(phi1)
            d1=norm(cross(or1-[phi1(n) psi(n) phi2(n)],or2-[phi1(n) psi(n) phi2(n)]))/norm(or1-or2);
            d2=norm(cross(or1-[phi1(n)-360 psi(n) phi2(n)],or2-[phi1(n)-360 psi(n) phi2(n)]))/norm(or1-or2);
            d3=norm(cross(or1-[phi1(n)-360 psi(n)-180 phi2(n)],or2-[phi1(n)-360 psi(n)-180 phi2(n)]))/norm(or1-or2);
            d4=norm(cross(or1-[phi1(n) psi(n)-180 phi2(n)],or2-[phi1(n) psi(n)-180 phi2(n)]))/norm(or1-or2);
            d5=norm(cross(or1-[phi1(n) psi(n) phi2(n)-360],or2-[phi1(n) psi(n) phi2(n)-180]))/norm(or1-or2);
            d6=norm(cross(or1-[phi1(n)-360 psi(n) phi2(n)-360],or2-[phi1(n)-360 psi(n) phi2(n)-180]))/norm(or1-or2);
            d7=norm(cross(or1-[phi1(n)-360 psi(n)-180 phi2(n)-360],or2-[phi1(n)-360 psi(n)-180 phi2(n)-180]))/norm(or1-or2);
            d8=norm(cross(or1-[phi1(n) psi(n)-180 phi2(n)-360],or2-[phi1(n) psi(n)-180 phi2(n)-180]))/norm(or1-or2);

            d(n)=min([d1 d2 d3 d4 d5 d6 d7 d8]);
        end

        odf=exp(-d.^2/(2*sigma)^2);
end

% odf=F/sum(F/(8*pi^2),'all'); % the normalisation is a question, one need to compute it