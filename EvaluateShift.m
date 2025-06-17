function [results,ahkl,ahkl_T,iserr]=EvaluateShift(data,GlobOpt,FitOpt)
hkl=[data.h data.k data.l];
sin2psi=(sind(data.psi)).^2;
tth=data.tth;
switch FitOpt.sys
    case {'cubic_fcc','cubic_bcc'}
        ahkl=GlobOpt.lam./(2*sind(tth/2)).*sqrt(hkl(:,1).^2+hkl(:,2).^2+hkl(:,3).^2);
        if ahkl~=-1
            [results,ahkl_T,iserr]=EvaluatePeakShift(sin2psi,ahkl,tth,hkl,GlobOpt,FitOpt);
        end
    case 'hexagonal'
        ahkl=GlobOpt.lam./(2*sind(tth/2)); % ahkl is dhkl, but the same notation is used in the scripts
        if ahkl~=-1
            [results,ahkl_T,iserr]=EvaluatePeakShift_hexagonal(sin2psi,ahkl,tth,hkl,GlobOpt,FitOpt);
        end
    otherwise
        disp('Calculation is not available for this system.');
        ahkl=-1;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% AUXILIARY FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=ahkl_shift(x,p,hkl,tth,GlobOpt,FitOpt)
% Kernel function for the ahkl shift fit.
% Depending on the model for the general stress influence, different number
% of fitting parameters are needed (p(1:4)). Therefore, some of them need
% to be fixed to 0 during the fit depending on the choosen model.
% The evaluation of a0 and other constants is done separately after fit,
% in dependence on the model and also in dependence on availability of
% material constants input by the user.
h=hkl(:,1);
k=hkl(:,2);
l=hkl(:,3);

% general stress influence
switch FitOpt.PeakShiftModel
    case 'NoStress'
        % p(1)=a0
        y_stress=p(1);
    case 'IsotropicStress'
        % a0*(sig*(1/2*s2*sin2psi+2*s1)+1)
        y_stress=p(1).*x+p(2);
    case 'VookWitt'
        Gamma=(h.^2.*k.^2+k.^2.*l.^2+h.^2.*l.^2)./((h.^2+k.^2+l.^2).^2);
        y_stress=p(1).*x+p(2).*Gamma.*x+p(3).*Gamma+p(4);
    case 'Reuss'
        Gamma=(h.^2.*k.^2+k.^2.*l.^2+h.^2.*l.^2)./((h.^2+k.^2+l.^2).^2);
        y_stress=p(1).*x-3*p(2).*Gamma.*x+2*p(2).*Gamma+p(3);
end

% stacking faults influence (powder)
% p(5)=a0*alpha
switch FitOpt.sys
    case 'cubic_fcc'
        y_SF=zeros(size(x));
        for n=1:numel(h)
            [planes,m_hkl]=GetEquivalentPlanes([h(n) k(n) l(n)],'cubic');
            L=planes(:,1)+planes(:,2)+planes(:,3);
            b=mod(L,3);
            b(b==2)=-1;
            y_SF(n)=p(5).*sqrt(3)/(4*pi).*sum(b.*L)./(m_hkl.*(h(n)^2+k(n)^2+l(n)^2));
        end
    case 'cubic_bcc'
        y_SF=zeros(size(x));
        % It is zero from definition.
        %         for n=1:numel(h)
        %             [planes,m_hkl]=GetEquivalentPlanes([h(n) k(n) l(n)],'cubic');
        %             L=-planes(:,1)-planes(:,2)+2*planes(:,3);
        %             b=mod(L,3);
        %             b(b==2)=-1;
        %             y_SF(n)=p(5).*sqrt(3)/(4*pi).*sum(b.*L)./(m_hkl.*sqrt(h(n)^2+k(n)^2+l(n)^2));
        %         end
    case 'hexagonal'

end

% vertical sample shift influence (divergent beam and Bragg-Brentano geometry)
% p(6)=a0*sycos
y_VSShift=p(6)./tand(tth/2);

% detector shift influence (divergent beam and Bragg-Brentano geometry)
% p(7)=a0*dth
y_DetShift=p(7).*cosd(tth/2)./tand(tth/2);

y=y_stress+y_SF+y_VSShift+y_DetShift;
end


function y=dhkl_shift_hexagonal(x,p,hkl,tth,GlobOpt,FitOpt)
% Kernel function for the ahkl shift fit.
% Depending on the model for the general stress influence, different number
% of fitting parameters are needed (p(1:4)). Therefore, some of them need
% to be fixed to 0 during the fit depending on the choosen model.
% The evaluation of a0 and other constants is done separately after fit,
% in dependence on the model and also in dependence on availability of
% material constants input by the user.
h=hkl(:,1);
k=hkl(:,2);
l=hkl(:,3);

% B=Reci(p(1),p(1),p(2),90,90,120);
% dhkl0=zeros(size(tth));
% for n=1:numel(tth)
%     dhkl0(n)=2*pi./norm(B*hkl(n,:)');
% end

% interplanar distance divided by a^2
dhkl0_=1./sqrt(4/3.*(h.^2+h.*k+k.^2)+p(2).^2.*l.^2);

% general stress influence
switch FitOpt.PeakShiftModel
    case 'NoStress'
        % p(1)=a0
        y_stress=p(1).*dhkl0_;
    case 'IsotropicStress'
        % a0*(sig*(1/2*s2*sin2psi+2*s1)+1)
        y_stress=dhkl0_.*(p(3).*x+p(4));
    case 'VookWitt'
%         Gamma=(h.^2.*k.^2+k.^2.*l.^2+h.^2.*l.^2)./((h.^2+k.^2+l.^2).^2);
%         y_stress=p(1).*x+p(2).*Gamma.*x+p(3).*Gamma+p(4);
    case 'Reuss'
%         Gamma=(h.^2.*k.^2+k.^2.*l.^2+h.^2.*l.^2)./((h.^2+k.^2+l.^2).^2);
%         y_stress=p(1).*x-3*p(2).*Gamma.*x+2*p(2).*Gamma+p(3);
end

% stacking faults influence (powder)
% p(5)=a0*alpha
% switch GlobOpt.sys
%     case 'hexagonal'
%         % there is no peak shift
% end

% vertical sample shift influence (divergent beam and Bragg-Brentano geometry)
% p(6)=a0*sycos
y_VSShift=p(6)./tand(tth/2);

% detector shift influence (divergent beam and Bragg-Brentano geometry)
% p(7)=a0*dth
y_DetShift=p(7).*cosd(tth/2)./tand(tth/2);

y=y_stress+y_VSShift+y_DetShift;
end


function y=ahkl_ResultsEvaluation(pT,ss,x,ahkl,hkl,GlobOpt,FitOpt)
h=hkl(:,1);
k=hkl(:,2);
l=hkl(:,3);
switch FitOpt.PeakShiftModel
    case 'NoStress'
        % p(1)=a0
        y.a0=pT(1);
        y.a0_err=ss(1);
    case 'IsotropicStress'
        % a0*(sig*(1/2*s2*sin2psi+2*s1)+1)
        y.a_norm=pT(2);
        y.a_norm_err=ss(2);
        y.a_par=pT(1)+pT(2);
        y.a_par_err=ss(1)+ss(2);
        y.sig_ef=pT(1); % a0*sig*(1+nu)/E
        y.sig_ef_err=ss(1); % a0*sig*(1+nu)/E
        if FitOpt.Material.av==1
            y.sin2psi0=2*FitOpt.Material.nu/(1+FitOpt.Material.nu);

            y.a0=(y.a_par-y.a_norm)*y.sin2psi0+y.a_norm;
            y.a0_err=(y.a_par_err-y.a_norm_err)*y.sin2psi0+y.a_norm_err;
            y.sig=y.sig_ef*FitOpt.Material.E/(y.a0.*(1+FitOpt.Material.nu));
            y.sig_err=y.sig_ef_err*FitOpt.Material.E/(y.a0.*(1+FitOpt.Material.nu))+y.sig_ef*y.a0_err*FitOpt.Material.E/(y.a0.^2.*(1+FitOpt.Material.nu));
        end
    case 'VookWitt'
        % first the anisotropy needs to be substracted
        y.a_norm=pT(4);
        y.a_norm_err=ss(4);
        y.a_par=pT(1)+pT(4);
        y.a_par_err=ss(1)+ss(4);
        y.sig_ef=pT(1); % a0*sig*(1+nu_h00)/E_h00, but the stress needs to be evaluated from p(1) from the original fit, for cubic crystals E_h00=E and nu_h00=nu
        y.sig_ef_err=ss(1); % a0*sig*(1+nu_h00)/E_h00, but the stress needs to be evaluated from p(1) from the original fit, for cubic crystals E_h00=E and nu_h00=nu
        if FitOpt.Material.av==1
            if isnan(FitOpt.Material.E_100)
                y.sin2psi0=2*FitOpt.Material.nu/(1+FitOpt.Material.nu);

                y.a0=(y.a_par-y.a_norm)*y.sin2psi0+y.a_norm;
                y.a0_err=(y.a_par_err-y.a_norm_err)*y.sin2psi0+y.a_norm_err;
                y.sig=y.sig_ef*FitOpt.Material.E/(y.a0.*(1+FitOpt.Material.nu));
                y.sig_err=y.sig_ef_err*FitOpt.Material.E/(y.a0.*(1+FitOpt.Material.nu))+y.sig_ef*y.a0_err*FitOpt.Material.E/(y.a0.^2.*(1+FitOpt.Material.nu));
            else
                y.sin2psi0=2*FitOpt.Material.nu_100/(1+FitOpt.Material.nu_100);

                y.a0=(y.a_par-y.a_norm)*y.sin2psi0+y.a_norm;
                y.a0_err=(y.a_par_err-y.a_norm_err)*y.sin2psi0+y.a_norm_err;
                y.sig=y.sig_ef*FitOpt.Material.E_100/(y.a0.*(1+FitOpt.Material.nu_100));
                y.sig_err=y.sig_ef_err*FitOpt.Material.E_100/(y.a0.*(1+FitOpt.Material.nu_100))+y.a0_err*FitOpt.Material.E_100/(y.a0.^2.*(1+FitOpt.Material.nu_100));
            end
        end
    case 'Reuss'
        % first the anisotropy needs to be substracted
        y.a_norm=pT(3);
        y.a_norm_err=ss(3);
        y.a_par=pT(1)+pT(3);
        y.a_par_err=ss(1)+ss(3);
        y.sig_ef=pT(1); % a0*sig*(1+nu_h00)/E_h00, but the stress needs to be evaluated from p(1) from the original fit, for cubic crystals E_h00=E and nu_h00=nu
        y.sig_ef_err=ss(1); % a0*sig*(1+nu_h00)/E_h00, but the stress needs to be evaluated from p(1) from the original fit, for cubic crystals E_h00=E and nu_h00=nu
        if FitOpt.Material.av==1
            if isnan(FitOpt.Material.E_100)
                y.sin2psi0=2*FitOpt.Material.nu/(1+FitOpt.Material.nu);

                y.a0=(y.a_par-y.a_norm)*y.sin2psi0+y.a_norm;
                y.a0_err=(y.a_par_err-y.a_norm_err)*y.sin2psi0+y.a_norm_err;
                y.sig=y.sig_ef*FitOpt.Material.E/(y.a0.*(1+FitOpt.Material.nu));
                y.sig_err=y.sig_ef_err*FitOpt.Material.E/(y.a0.*(1+FitOpt.Material.nu))+y.sig_ef*y.a0_err*FitOpt.Material.E/(y.a0.^2.*(1+FitOpt.Material.nu));
            else
                y.sin2psi0=2*FitOpt.Material.nu_100/(1+FitOpt.Material.nu_100);

                y.a0=(y.a_par-y.a_norm)*y.sin2psi0+y.a_norm;
                y.a0_err=(y.a_par_err-y.a_norm_err)*y.sin2psi0+y.a_norm_err;
                y.sig=y.sig_ef*FitOpt.Material.E_100/(y.a0.*(1+FitOpt.Material.nu_100));
                y.sig_err=y.sig_ef_err*FitOpt.Material.E_100/(y.a0.*(1+FitOpt.Material.nu_100))+y.a0_err*FitOpt.Material.E_100/(y.a0.^2.*(1+FitOpt.Material.nu_100));
            end
        end
end

if FitOpt.FitSF==1
    % p(5)=a0*alpha
    y.alphaSF_ef=pT(5);
    if FitOpt.Material.av==1
       y.alphaSF=pT(5)./y.a0;
    end
end

if GlobOpt.FitVerticalShift==1
    % p(6)=a0*sycos=a0*s/2R
    y.sycos_ef=pT(6);
    y.sycos_ef_err=ss(6);
    if FitOpt.Material.av==1
       y.sycos=pT(6)./y.a0;
       y.sycos_err=ss(6)./y.a0;
    end
end

if GlobOpt.FitDetectorShift==1
    % p(7)=a0*dth
    y.dth_ef=pT(7);
    y.dth_ef_err=ss(7);
    if FitOpt.Material.av==1
       y.dth=pT(7)./y.a0;
       y.dth_err=ss(7)./y.a0;
    end
end
end


function y=dhkl_ResultsEvaluation_hexagonal(pT,ss,x,ahkl,hkl,GlobOpt,FitOpt)
h=hkl(:,1);
k=hkl(:,2);
l=hkl(:,3);
switch FitOpt.PeakShiftModel
    case 'NoStress'
        % p(1)=a0, p(2)=a0/c0
        y.a0=pT(1);
        y.a0_err=ss(1);
        y.acrat=pT(2);
        y.acrat_err=ss(2);
        y.c0=pT(1)/pT(2);
        y.c0_err=sqrt((ss(1)/pT(2))^2+(pT(1)*ss(2)/(pT(2))^2)^2);
    case 'IsotropicStress'
        y.a_norm=pT(4);
        y.a_norm_err=ss(4);
        y.a_par=pT(3)+pT(4);
        y.a_par_err=ss(3)+ss(4);
        % a0*(sig*(1/2*s2*sin2psi+2*s1)+1)
        y.acrat=pT(2);
        y.acrat_err=ss(2);
        y.sig_ef=pT(3); % a0*sig*(s11+s12-2*s13)
        y.sig_ef_err=ss(3); % a0*sig*(s11+s12-2*s13)
        y.s13term=pT(4); % 2*a0*sig*s13+1
        y.s13term_err=ss(4); % 2*a0*sig*s13+1
        if FitOpt.Material.av==1
            y.sin2psi0=2*FitOpt.Material.nu/(1+FitOpt.Material.nu);

            y.a0=(y.a_par-y.a_norm)*y.sin2psi0+y.a_norm;
            y.a0_err=(y.a_par_err-y.a_norm_err)*y.sin2psi0+y.a_norm_err;
            y.sig=y.sig_ef*FitOpt.Material.E/(y.a0.*(1+FitOpt.Material.nu));
            y.sig_err=y.sig_ef_err*FitOpt.Material.E/(y.a0.*(1+FitOpt.Material.nu))+y.sig_ef*y.a0_err*FitOpt.Material.E/(y.a0.^2.*(1+FitOpt.Material.nu));
            y.c0=y.a0/y.acrat;
            y.c0_err=sqrt((y.a0_err/y.acrat)^2+(y.a0*y.acrat_err/(y.acrat)^2)^2);
        end
    case 'VookWitt'
%         % first the anisotropy needs to be substracted
%         y.a_norm=pT(4);
%         y.a_norm_err=ss(4);
%         y.a_par=pT(1)+pT(4);
%         y.a_par_err=ss(1)+ss(4);
%         y.sig_ef=pT(1); % a0*sig*(1+nu_h00)/E_h00, but the stress needs to be evaluated from p(1) from the original fit, for cubic crystals E_h00=E and nu_h00=nu
%         y.sig_ef_err=ss(1); % a0*sig*(1+nu_h00)/E_h00, but the stress needs to be evaluated from p(1) from the original fit, for cubic crystals E_h00=E and nu_h00=nu
%         if GlobOpt.Material.av==1
%             if isnan(GlobOpt.Material.E_100)
%                 y.sin2psi0=2*GlobOpt.Material.nu/(1+GlobOpt.Material.nu);
% 
%                 y.a0=(y.a_par-y.a_norm)*y.sin2psi0+y.a_norm;
%                 y.a0_err=(y.a_par_err-y.a_norm_err)*y.sin2psi0+y.a_norm_err;
%                 y.sig=y.sig_ef*GlobOpt.Material.E/(y.a0.*(1+GlobOpt.Material.nu));
%                 y.sig_err=y.sig_ef_err*GlobOpt.Material.E/(y.a0.*(1+GlobOpt.Material.nu))+y.sig_ef*y.a0_err*GlobOpt.Material.E/(y.a0.^2.*(1+GlobOpt.Material.nu));
%             else
%                 y.sin2psi0=2*GlobOpt.Material.nu_100/(1+GlobOpt.Material.nu_100);
% 
%                 y.a0=(y.a_par-y.a_norm)*y.sin2psi0+y.a_norm;
%                 y.a0_err=(y.a_par_err-y.a_norm_err)*y.sin2psi0+y.a_norm_err;
%                 y.sig=y.sig_ef*GlobOpt.Material.E_100/(y.a0.*(1+GlobOpt.Material.nu_100));
%                 y.sig_err=y.sig_ef_err*GlobOpt.Material.E_100/(y.a0.*(1+GlobOpt.Material.nu_100))+y.a0_err*GlobOpt.Material.E_100/(y.a0.^2.*(1+GlobOpt.Material.nu_100));
%             end
%         end
    case 'Reuss'
%         % first the anisotropy needs to be substracted
%         y.a_norm=pT(3);
%         y.a_norm_err=ss(3);
%         y.a_par=pT(1)+pT(3);
%         y.a_par_err=ss(1)+ss(3);
%         y.sig_ef=pT(1); % a0*sig*(1+nu_h00)/E_h00, but the stress needs to be evaluated from p(1) from the original fit, for cubic crystals E_h00=E and nu_h00=nu
%         y.sig_ef_err=ss(1); % a0*sig*(1+nu_h00)/E_h00, but the stress needs to be evaluated from p(1) from the original fit, for cubic crystals E_h00=E and nu_h00=nu
%         if GlobOpt.Material.av==1
%             if isnan(GlobOpt.Material.E_100)
%                 y.sin2psi0=2*GlobOpt.Material.nu/(1+GlobOpt.Material.nu);
% 
%                 y.a0=(y.a_par-y.a_norm)*y.sin2psi0+y.a_norm;
%                 y.a0_err=(y.a_par_err-y.a_norm_err)*y.sin2psi0+y.a_norm_err;
%                 y.sig=y.sig_ef*GlobOpt.Material.E/(y.a0.*(1+GlobOpt.Material.nu));
%                 y.sig_err=y.sig_ef_err*GlobOpt.Material.E/(y.a0.*(1+GlobOpt.Material.nu))+y.sig_ef*y.a0_err*GlobOpt.Material.E/(y.a0.^2.*(1+GlobOpt.Material.nu));
%             else
%                 y.sin2psi0=2*GlobOpt.Material.nu_100/(1+GlobOpt.Material.nu_100);
% 
%                 y.a0=(y.a_par-y.a_norm)*y.sin2psi0+y.a_norm;
%                 y.a0_err=(y.a_par_err-y.a_norm_err)*y.sin2psi0+y.a_norm_err;
%                 y.sig=y.sig_ef*GlobOpt.Material.E_100/(y.a0.*(1+GlobOpt.Material.nu_100));
%                 y.sig_err=y.sig_ef_err*GlobOpt.Material.E_100/(y.a0.*(1+GlobOpt.Material.nu_100))+y.a0_err*GlobOpt.Material.E_100/(y.a0.^2.*(1+GlobOpt.Material.nu_100));
%             end
%         end
end

% if GlobOpt.FitSF==1
%     % p(5)=a0*alpha
%     y.alphaSF_ef=pT(5);
%     if GlobOpt.Material.av==1
%        y.alphaSF=pT(5)./y.a0;
%     end
% end

if GlobOpt.FitVerticalShift==1
    % p(5)=a0*sycos=a0*s/2R
    y.sycos_ef=pT(5);
    y.sycos_ef_err=ss(5);
    if FitOpt.Material.av==1
       y.sycos=pT(5)./y.a0;
       y.sycos_err=ss(5)./y.a0;
    end
end

if GlobOpt.FitDetectorShift==1
    % p(6)=a0*dth
    y.dth_ef=pT(6);
    y.dth_ef_err=ss(6);
    if FitOpt.Material.av==1
       y.dth=pT(6)./y.a0;
       y.dth_err=ss(6)./y.a0;
    end
end
end


function [results,ahkl_T,iserr]=EvaluatePeakShift(sin2psi,ahkl,tth,hkl,GlobOpt,FitOpt)
p0=zeros(1,7);
fix=ones(1,7);
lb=-Inf*ones(1,7);
ub=+Inf*ones(1,7);
switch FitOpt.PeakShiftModel
    case 'NoStress'
        p0(1)=mean(ahkl);
        fix(1:4)=[0 1 1 1];
        lb(1)=0.9*p0(1);
        ub(1)=1.1*p0(1);
    case 'IsotropicStress'
        pT0=polyfit(sin2psi,ahkl,1);
        p0(1:2)=[pT0(2) pT0(1)];
        fix(1:4)=[0 0 1 1];
%         lb(1:2)=[0.8*p0(1) -abs(p0(2))*50];
%         ub(1:2)=[1.2*p0(1) abs(p0(2))*50];
    case 'VookWitt'
        pT0=polyfit(sin2psi,ahkl,1);
        p0(1:2)=pT0;
        fix(1:4)=[0 0 0 0];
    case 'Reuss'
        pT0=polyfit(sin2psi,ahkl,1);
        p0(1:2)=pT0;
        fix(1:4)=[0 0 0 1];
end

if FitOpt.FitSF==1
    if strcmp(FitOpt.sys,'cubic_fcc') % for bcc and hcp there is no influence on peak shift
        p0(5)=1e-2; % it is zero by default
        fix(5)=0;
        lb(5)=0;
        ub(5)=1;
    end
end

if GlobOpt.FitVerticalShift==1
    p0(6)=0;
    fix(6)=0;
    lb(6)=-0.1; % in radians
    ub(6)=0.1; % in radians
end

if GlobOpt.FitDetectorShift==1
    p0(7)=0;
    fix(7)=0;
    lb(7)=-0.1; % in radians
    ub(7)=0.1; % in radians
end

if numel(tth)<=sum(fix==0)
    iserr=1; % more parameters than datapoints
elseif sum(imag(p0))~=0
    iserr=2; % initial guess gives complex parameters
else
    iserr=0;
    [pT1,resnorm,~,~,~,~,jacob]=lsqcurvefit(@(p,xdata) ahkl_shift(xdata,interlace(p0,p,fix),hkl,tth,GlobOpt,FitOpt),p0(~fix),sin2psi,ahkl,lb(~fix),ub(~fix));
    ss1=full(diag(sqrt(resnorm*inv(jacob'*jacob))/(length(sin2psi)-length(pT1))));
    pT=interlace(p0,pT1,fix); % for the correct evaluation it is necessary to have all parameters, included those fixed to zero
    ss=interlace(zeros(size(pT)),ss1,fix);

    if sum(imag([pT1(:); ss1(:)]))~=0 || max(abs(ss1(:)./pT1(:)))>0.5
        iserr=-1; % complex values of fitting parameters, or large errors
    end
end

if iserr<=0
    results=ahkl_ResultsEvaluation(pT,ss,sin2psi,ahkl,hkl,GlobOpt,FitOpt);
    ahkl_T=ahkl_shift(sin2psi,pT,hkl,tth,GlobOpt,FitOpt);
else
    results=NaN;
    ahkl_T=NaN;
end

end


function [results,dhkl_T,iserr]=EvaluatePeakShift_hexagonal(sin2psi,dhkl,tth,hkl,GlobOpt,FitOpt)
p0=zeros(1,7);
fix=ones(1,7);
lb=-Inf*ones(1,7);
ub=+Inf*ones(1,7);
[a0,c0]=EstimateHexagonalLattice(hkl,dhkl);
switch FitOpt.PeakShiftModel
    case 'NoStress'
        p0(1:2)=[a0 a0/c0];
        fix(1:4)=[0 0 1 1];
        lb(1:2)=[0.8*p0(1) 0];
        ub(1:2)=[1.2*p0(1) p0(2)*10];
    case 'IsotropicStress'
        p0(1:2)=[a0 a0/c0];
        lb(1:2)=[0.8*p0(1) 0];
        ub(1:2)=[1.2*p0(1) p0(2)*10];
        if FitOptOpt.FitHcpLP==1
            fix(1:4)=[1 0 0 0];
        else
            fix(1:4)=[1 1 0 0];
        end
    case 'VookWitt'
        disp('Sorry, I can''t handle this method yet.')
%         pT0=polyfit(sin2psi,ahkl,1);
%         p0(1:2)=pT0;
%         fix(1:4)=[0 0 0 0];
    case 'Reuss'
        disp('Sorry, I can''t handle this method yet.')
%         pT0=polyfit(sin2psi,ahkl,1);
%         p0(1:2)=pT0;
%         fix(1:4)=[0 0 0 1];
end

% % There is no peak shift
% if GlobOpt.FitSF==1
%     p0(6)=1e-2; % it is zero by default
%     fix(6)=0;
%     lb(6)=0;
%     ub(6)=1; % is the alpha density or probability?
% end

if GlobOpt.FitVerticalShift==1
    p0(5)=0;
    fix(5)=0;
    lb(5)=-0.1; % in radians
    ub(5)=0.1; % in radians
end

if GlobOpt.FitDetectorShift==1
    p0(6)=0;
    fix(6)=0;
    lb(6)=-0.1; % in radians
    ub(6)=0.1; % in radians
end

if numel(tth)<=sum(fix==0)
    iserr=1; % more parameters than datapoints
elseif sum(imag(p0))~=0
    iserr=2; % initial guess gives complex parameters
else
    iserr=0;
    [pT1,resnorm,~,~,~,~,jacob]=lsqcurvefit(@(p,xdata) dhkl_shift_hexagonal(xdata,interlace(p0,p,fix),hkl,tth,GlobOpt,FitOpt),p0(~fix),sin2psi,dhkl,lb(~fix),ub(~fix));
    ss1=full(diag(sqrt(resnorm*inv(jacob'*jacob))/(length(sin2psi)-length(pT1))));
    pT=interlace(p0,pT1,fix); % for the correct evaluation it is necessary to have all parameters, included those fixed to zero
    ss=interlace(zeros(size(pT)),ss1,fix);

    if sum(imag([pT1(:); ss1(:)]))~=0 || max(abs(ss1(:)./pT1(:)))>0.5
        iserr=-1; % complex values of fitting parameters, or large errors
    end
end

if iserr<=0
    results=dhkl_ResultsEvaluation_hexagonal(pT,ss,sin2psi,dhkl,hkl,GlobOpt,FitOpt);
    dhkl_T=dhkl_shift_hexagonal(sin2psi,pT,hkl,tth,GlobOpt,FitOpt);
else
    results=NaN;
    dhkl_T=NaN;
end

end