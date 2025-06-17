function [E_hkl,nu_hkl,G_hkl]=ElasticConstantsHKL(Sij,h,k,l)
% h=hkl(:,1);
% k=hkl(:,2);
% l=hkl(:,3);
% switch sys
%     case 'cubic'
        Gamma=(h.^2.*k.^2+h.^2.*l.^2+k.^2.*l.^2)./(h.^2+k.^2+l.^2).^2;
        sA=Sij(1,1)-Sij(1,2)-0.5*Sij(4,4);
        E_hkl=1./(Sij(1,1)-2*sA.*Gamma);
        G_hkl=1./(Sij(4,4)+4*((Sij(1,1)-Sij(1,2))-0.5*Sij(4,4)).*Gamma);
%         nu_hkl=E_hkl./(2*G_hkl)-1;

%         sA=Sij(1,1)-Sij(1,2)-0.5*Sij(4,4);
        ind=find(h==k & k==0);
        sp13=Sij(1,2)+sA.*(1+(h.^4+k.^4)./((h.^2+k.^2).^2)).*((h.*l).^2+(k.*l).^2)./((h.^2+k.^2+l.^2).^2);
        sp23=Sij(1,2)+2*sA.*(h.^2.*k.^2)./((h.^2+k.^2).^2).*(h.^2+k.^2)./(h.^2+k.^2+l.^2);
        sp13(ind)=Sij(1,2);%+sA.*(1+1).*(h(ind).^2.*l(ind).^2+k(ind).^2.*l(ind).^2)./(h(ind).^2+k(ind).^2+l(ind).^2).^2;
        sp23(ind)=Sij(1,2);%+2*sA.*0.^2.*(h(ind).^2+k(ind).^2)./(h(ind).^2+k(ind).^2+l(ind).^2).^2;

%         cal=h./sqrt(h.^2+k.^2+l.^2);
%         cbe=k./sqrt(h.^2+k.^2+l.^2);
%         cga=l./sqrt(h.^2+k.^2+l.^2);
%         psi=acosd(cga);
%         phi=atan2d(cbe,cal);
%         nu_hkl=zeros(size(h));
%         for n=1:numel(h)
%             Sij_rot=GetReducedComplianceTensor(RotateFullElasticTensor(GetFullComplianceTensor(Sij),RotMat(phi(n),psi(n),0,'zxz')));
%             nu_hkl(n)=-0.5*(Sij_rot(1,2)+Sij_rot(1,3))/Sij_rot(1,1)/3+...
%                 -0.5*(Sij_rot(2,1)+Sij_rot(2,3))/Sij_rot(2,2)/3+...
%                 -0.5*(Sij_rot(3,1)+Sij_rot(3,2))/Sij_rot(3,3)/3;
%         end
% 
%         nu_hkl=(Sij(1,2)+sA./(h.^2+k.^2+l.^2).*((h.^2.*l./(sqrt(h.^2+k.^2).*sqrt(h.^2+k.^2+l.^2))).^2+(k.^2.*l./(sqrt(h.^2+k.^2).*sqrt(h.^2+k.^2+l.^2))).^2+(l.*sqrt(h.^2+k.^2)./sqrt(h.^2+k.^2+l.^2)).^2))./(-Sij(1,1)+2*sA.*Gamma);
%         nu_hkl(ind)=-2*Sij(1,2)./(2*Sij(1,1)-4*sA.*Gamma(ind));
%         nu_hkl=-(2*Sij(1,2)+sA.*((1+(h.^4+k.^4)./(h.^2+k.^2).^2).*l.^2./(h.^2+k.^2+l.^2)+2*h.^2.*k.^2./(h.^2+k.^2).^2).*(h.^2+k.^2)./(h.^2+k.^2+l.^2))./(2*Sij(1,1)-4*sA.*Gamma);
        nu_hkl=-(sp13+sp23)./(2*Sij(1,1)-4*sA.*Gamma);

%         s1=-nu_hkl./E_hkl;
%         s2=2*(1+nu_hkl)./E_hkl;
%     case 'hexagonal'
%         ca3=sqrt(hkl(:,3).^2./(hkl(:,1).^2+hkl(:,2).^2+hkl(:,3).^2).^2);
%         E_hkl=1./(Sij(1,1)*(1-ca3.^2).^2+Sij(3,3)*ca3.^4+(2*Sij(1,3)+Sij(4,4)).*ca3.^2.*(1-ca3.^2));
%         G_hkl=1./(Sij(4,4)+(Sij(1,1)-Sij(1,2)-0.5*Sij(4,4))*(1-ca3.^2)+2*(Sij(1,1)+Sij(3,3)-2*Sij(1,3)-Sij(4,4)).*ca3.^2.*(1-ca3.^2));
%         nu_hkl=E_hkl./(2*G_hkl)-1;
% 
%         s1=-nu_hkl./E_hkl;
%         s2=2*(1+nu_hkl)./E_hkl;
%     otherwise
%         s1=NaN;
%         s2=NaN;
%         errordlg('XEC model for this system is unfortunately not available.');
% end