function [s1,s2]=XrayElasticConstants(hkl,Sij,sys,XECmodel,KroenerCoeff)
h=hkl(:,1);
k=hkl(:,2);
l=hkl(:,3);
switch sys
    case 'cubic'
        Gamma=(h.^2.*k.^2+h.^2.*l.^2+k.^2.*l.^2)./(h.^2+k.^2+l.^2).^2;
        sA=Sij(1,1)-Sij(1,2)-0.5*Sij(4,4);
        E_hkl=1./(Sij(1,1)-2*sA.*Gamma);
%         G_hkl=1./(Sij(4,4)+4*((Sij(1,1)-Sij(1,2))-0.5*Sij(4,4)).*Gamma);
%         nu_hkl=E_hkl./(2*G_hkl)-1;

        sA=Sij(1,1)-Sij(1,2)-0.5*Sij(4,4);
        ind=find(h==k & k==0);
        sp13=Sij(1,2)+sA.*(1+(h.^4+k.^4)./((h.^2+k.^2).^2)).*((h.*l).^2+(k.*l).^2)./((h.^2+k.^2+l.^2).^2);
        sp23=Sij(1,2)+2*sA.*(h.^2.*k.^2)./((h.^2+k.^2).^2).*(h.^2+k.^2)./(h.^2+k.^2+l.^2);
        sp13(ind)=Sij(1,2);%+sA.*(1+1).*(h(ind).^2.*l(ind).^2+k(ind).^2.*l(ind).^2)./(h(ind).^2+k(ind).^2+l(ind).^2).^2;
        sp23(ind)=Sij(1,2);%+2*sA.*0.^2.*(h(ind).^2+k(ind).^2)./(h(ind).^2+k(ind).^2+l(ind).^2).^2;

%         nu_hkl=(Sij(1,2)+sA./(h.^2+k.^2+l.^2).*((h.^2.*l./(sqrt(h.^2+k.^2).*sqrt(h.^2+k.^2+l.^2))).^2+(k.^2.*l./(sqrt(h.^2+k.^2).*sqrt(h.^2+k.^2+l.^2))).^2+(l.*sqrt(h.^2+k.^2)./sqrt(h.^2+k.^2+l.^2)).^2))./(-Sij(1,1)+2*sA.*Gamma);
%         nu_hkl(ind)=-2*Sij(1,2)./(2*Sij(1,1)-4*sA.*Gamma(ind));
%         nu_hkl=-(2*Sij(1,2)+sA.*((1+(hkl(:,1).^4+hkl(:,2).^4)./(hkl(:,1).^2+hkl(:,2).^2).^2).*hkl(:,3).^2./(hkl(:,1).^2+hkl(:,2).^2+hkl(:,3).^2)+2*hkl(:,1).^2.*hkl(:,2).^2./(hkl(:,1).^2+hkl(:,2).^2).^2).*(hkl(:,1).^2+hkl(:,2).^2)./(hkl(:,1).^2+hkl(:,2).^2+hkl(:,3).^2))./(2*Sij(1,1)-4*sA.*Gamma);
        nu_hkl=-(sp13+sp23)./(2*Sij(1,1)-4*sA.*Gamma);

        s1=-nu_hkl./E_hkl;
        s2=2*(1+nu_hkl)./E_hkl;

        s1_Voigt=(2*sA*(Sij(1,1)+2*Sij(1,2))+5*Sij(1,2)*Sij(4,4))/(6*sA+5*Sij(4,4))*ones(numel(h),1);
        s2_Voigt=2*(5*(Sij(1,1)-Sij(1,2))*Sij(4,4))/(6*sA+5*Sij(4,4))*ones(numel(h),1);
        s1_Reuss=Sij(1,2)+sA.*Gamma;
        s2_Reuss=2*(Sij(1,1)-Sij(1,2)-3*sA.*Gamma);
    case 'hexagonal'
        ca3=sqrt(hkl(:,3).^2./(hkl(:,1).^2+hkl(:,2).^2+hkl(:,3).^2).^2);
        E_hkl=1./(Sij(1,1)*(1-ca3.^2).^2+Sij(3,3)*ca3.^4+(2*Sij(1,3)+Sij(4,4)).*ca3.^2.*(1-ca3.^2));
        G_hkl=1./(Sij(4,4)+(Sij(1,1)-Sij(1,2)-0.5*Sij(4,4))*(1-ca3.^2)+2*(Sij(1,1)+Sij(3,3)-2*Sij(1,3)-Sij(4,4)).*ca3.^2.*(1-ca3.^2));
        nu_hkl=E_hkl./(2*G_hkl)-1;

        s1=-nu_hkl./E_hkl;
        s2=2*(1+nu_hkl)./E_hkl;

        s1_Voigt=(2*sA*(Sij(1,1)+2*Sij(1,2))+5*Sij(1,2)*Sij(4,4))/(6*sA+5*Sij(4,4))*ones(numel(h),1);
        s2_Voigt=2*(5*(Sij(1,1)-Sij(1,2))*Sij(4,4))/(6*sA+5*Sij(4,4))*ones(numel(h),1);
        s1_Reuss=-nu_hkl./E_hkl;
        s2_Reuss=2*(1+nu_hkl)./E_hkl;
    otherwise
        s1_Voigt=NaN;
        s2_Voigt=NaN;
        s1_Reuss=NaN;
        s2_Reuss=NaN;
        errordlg('XEC model for this system is unfortunately not available.');
end

switch XECmodel
    case 'Voigt'
        s1=s1_Voigt;
        s2=s2_Voigt;
    case 'Reuss'
        s1=s1_Reuss;
        s2=s2_Reuss;
    case 'Kroener'
        s1=s1_Voigt*(1-KroenerCoeff)+s1_Reuss*KroenerCoeff;
        s2=s2_Voigt*(1-KroenerCoeff)+s2_Reuss*KroenerCoeff;
        KroenerCoeff
        disp([s1_Voigt s1_Reuss s1])
        disp([s2_Voigt s2_Reuss s2])
    otherwise
        s1=NaN;
        s2=NaN;
        errordlg('Uknown model for XECs.')
end