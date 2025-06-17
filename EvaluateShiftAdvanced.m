function [results,ahkl,ahkl_T,iserr]=EvaluateShiftAdvanced(data,GlobOpt,FitOpt,Sij,Cij,ng)
hkl=[data.h data.k data.l];
sin2psi=(sind(data.psi)).^2;
tth=data.tth;
psi=data.psi;
phi=data.phi;
AdvancedShiftOptions=FitOpt.AdvancedShiftOptions;

switch FitOpt.sys
    case {'cubic_fcc','cubic_bcc'}
        ahkl=GlobOpt.lam./(2*sind(tth/2)).*sqrt(hkl(:,1).^2+hkl(:,2).^2+hkl(:,3).^2);
        disp(ahkl)
        if ahkl~=-1
            if strcmp(AdvancedShiftOptions.MaterialInfo,'StressFactors') 
                if strcmp(AdvancedShiftOptions.ODFtype,'PolycrystalWithTexture')
                    F=StressFactorsEvaluation(Sij,Cij,AdvancedShiftOptions.ODFfilename,AdvancedShiftOptions.StressFactorsMethod,hkl,psi,phi,AdvancedShiftOptions.Geometry,'cubic',NaN,0);
                elseif strcmp(AdvancedShiftOptions.ODFtype,'IdealPolycrystal')
                    F=StressFactorsEvaluation(Sij,Cij,1,AdvancedShiftOptions.StressFactorsMethod,hkl,psi,phi,AdvancedShiftOptions.Geometry,'cubic',NaN,0);
                end
            elseif strcmp(AdvancedShiftOptions.MaterialInfo,'XEC') 
                [s1,s2]=XrayElasticConstants(hkl,Sij,'cubic',AdvancedShiftOptions.XECmodel,AdvancedShiftOptions.KroenerCoeff);
                F=[s1,s2];
            else
                % in the cases of single-crystalline variants, F could be
                % computed relatively easily

                % compute rotation matrix
                EA=[AdvancedShiftOptions.OrientationAngles{ng,1} AdvancedShiftOptions.OrientationAngles{ng,2} AdvancedShiftOptions.OrientationAngles{ng,3}];
                rotopt1=[AdvancedShiftOptions.OrientationAngles{ng,4} AdvancedShiftOptions.OrientationAngles{ng,5} AdvancedShiftOptions.OrientationAngles{ng,6}];
                rotopt2=AdvancedShiftOptions.OrientationAngles{ng,7};
                R=RotMat(EA(1),EA(2),EA(3),rotopt1,rotopt2);

                % rotate compliance tensor
                Sfull=GetFullComplianceTensor(Sij);
                Sfull_rot=RotateFullElasticTensor(Sfull,R);

                % sum over 'measurement projections'
                F=zeros(3,3,numel(psi));
                for m=1:numel(psi)
                    switch AdvancedShiftOptions.Geometry
                        case 'GAXRD'
                            measd=[sind(psi(m)).*cosd(phi(m)); sind(psi(m)).*sind(phi(m)); cosd(psi(m))];
                        case 'Pole figure'
                            measd=[sind(psi(m)).*sind(phi(m)); -sind(psi(m)).*cosd(phi(m)); cosd(psi(m))];
%                             measd=[sind(psi(m)); sind(psi(m)); cosd(psi(m))];
                    end
                    
                    % compute stress factors from the rotated compliance
                    % matrix/tensor
                    for i=1:3
                        for j=1:3
                            for k=1:3
                                for l=1:3
                                    F(i,j,m)=F(i,j,m)+Sfull_rot(i,j,k,l)*measd(k)*measd(l);
%                                     F(i,j,m)=F(i,j,m)+Sfull_rot(i,j,k,l)*measd(k)*measd(l);
                                end
                            end
                        end
                    end
                end
            end
%             disp(F)
%             disp(ahkl)
%             disp(sin2psi)
%             disp(hkl)
%             disp(tth)
%             disp(psi)
%             disp(phi)
            [results,ahkl_T,iserr]=EvaluatePeakShiftAdvanced(F,ahkl,sin2psi,hkl,tth,psi,phi,GlobOpt,FitOpt,AdvancedShiftOptions);
        end
    case 'hexagonal'
        ahkl=GlobOpt.lam./(2*sind(tth/2)); % ahkl is dhkl, but the same notation is used in the scripts
        if ahkl~=-1
            [results,ahkl_T,iserr]=EvaluatePeakShiftAdvanced_hexagonal(sin2psi,ahkl,hkl,tth,psi,phi,GlobOpt,FitOpt,AdvancedShiftOptions);
        end
    otherwise
        disp('Calculation is not available for this system.');
        ahkl=-1;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% AUXILIARY FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=ahkl_shift_advanced(x,par,F,hkl,tth,psi,phi,GlobOpt,FitOpt,AdvancedShiftOptions)
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
sigma=zeros(6,1);
p=par(5:4+AdvancedShiftOptions.FitN); % neccessary to make the following eval work
for n=1:6
    sigma(n)=eval([AdvancedShiftOptions.StressComponents{n} ';']);
end
sig_ij=[sigma(1) sigma(6) sigma(5);
        sigma(6) sigma(2) sigma(4);
        sigma(5) sigma(4) sigma(3);];

eps=zeros(numel(tth),1);
switch AdvancedShiftOptions.MaterialInfo
    case 'XEC'
        eps=0.5*F(:,2).*(sig_ij(1,1).*cosd(phi).^2+sig_ij(1,2).*sind(2*phi)+sig_ij(2,2).*sind(phi).^2).*sind(psi).^2+...
            0.5*F(:,2).*(sig_ij(1,3).*cosd(phi)+sig_ij(2,3).*sind(phi)).*sind(2*psi)+...
            0.5*F(:,2).*sig_ij(3,3).*cosd(psi).^2+...
            F(:,1).*(sig_ij(1,1)+sig_ij(2,2)+sig_ij(3,3));
    otherwise
        for m=1:numel(tth)
            hlp=F(:,:,m).*sig_ij;
            eps(m)=sum(hlp,'all');
        end
end

y_stress=par(1).*(eps+1);

% stacking faults influence (powder)
% par(2)=a0*alpha
switch FitOpt.sys
    case 'cubic_fcc'
        y_SF=zeros(size(x));
        for n=1:numel(h)
            [planes,m_hkl]=GetEquivalentPlanes([h(n) k(n) l(n)],'cubic');
            L=planes(:,1)+planes(:,2)+planes(:,3);
            b=mod(L,3);
            b(b==2)=-1;
            y_SF(n)=par(2).*sqrt(3)/(4*pi).*sum(b.*L)./(m_hkl.*(h(n)^2+k(n)^2+l(n)^2));
        end
    case 'cubic_bcc'
        y_SF=zeros(size(x));
        % It is zero from definition.
        %         for n=1:numel(h)
        %             [planes,m_hkl]=GetEquivalentPlanes([h(n) k(n) l(n)],'cubic');
        %             L=-planes(:,1)-planes(:,2)+2*planes(:,3);
        %             b=mod(L,3);
        %             b(b==2)=-1;
        %             y_SF(n)=par(2).*sqrt(3)/(4*pi).*sum(b.*L)./(m_hkl.*sqrt(h(n)^2+k(n)^2+l(n)^2));
        %         end
    case 'hexagonal'

end

% vertical sample shift influence (divergent beam and Bragg-Brentano geometry)
% par(3)=a0*sycos
y_VSShift=par(3)./tand(tth/2);

% detector shift influence (divergent beam and Bragg-Brentano geometry)
% par(4)=a0*dth
y_DetShift=par(4).*cosd(tth/2)./tand(tth/2);

y=y_stress+y_SF+y_VSShift+y_DetShift;
end


% function y=dhkl_shift_hexagonal(x,p,hkl,tth,GlobOpt)
% % Kernel function for the ahkl shift fit.
% % Depending on the model for the general stress influence, different number
% % of fitting parameters are needed (p(1:4)). Therefore, some of them need
% % to be fixed to 0 during the fit depending on the choosen model.
% % The evaluation of a0 and other constants is done separately after fit,
% % in dependence on the model and also in dependence on availability of
% % material constants input by the user.
% h=hkl(:,1);
% k=hkl(:,2);
% l=hkl(:,3);
% 
% % B=Reci(p(1),p(1),p(2),90,90,120);
% % dhkl0=zeros(size(tth));
% % for n=1:numel(tth)
% %     dhkl0(n)=2*pi./norm(B*hkl(n,:)');
% % end
% 
% % interplanar distance divided by a^2
% dhkl0_=1./sqrt(4/3.*(h.^2+h.*k+k.^2)+p(2).^2.*l.^2);
% 
% % general stress influence
% switch GlobOpt.PeakShiftModel
%     case 'NoStress'
%         % p(1)=a0
%         y_stress=p(1).*dhkl0_;
%     case 'IsotropicStress'
%         % a0*(sig*(1/2*s2*sin2psi+2*s1)+1)
%         y_stress=dhkl0_.*(p(3).*x+p(4));
%     case 'VookWitt'
% %         Gamma=(h.^2.*k.^2+k.^2.*l.^2+h.^2.*l.^2)./((h.^2+k.^2+l.^2).^2);
% %         y_stress=p(1).*x+p(2).*Gamma.*x+p(3).*Gamma+p(4);
%     case 'Reuss'
% %         Gamma=(h.^2.*k.^2+k.^2.*l.^2+h.^2.*l.^2)./((h.^2+k.^2+l.^2).^2);
% %         y_stress=p(1).*x-3*p(2).*Gamma.*x+2*p(2).*Gamma+p(3);
% end
% 
% % stacking faults influence (powder)
% % p(5)=a0*alpha
% % switch GlobOpt.sys
% %     case 'hexagonal'
% %         % there is no peak shift
% % end
% 
% % vertical sample shift influence (divergent beam and Bragg-Brentano geometry)
% % p(6)=a0*sycos
% y_VSShift=p(6)./tand(tth/2);
% 
% % detector shift influence (divergent beam and Bragg-Brentano geometry)
% % p(7)=a0*dth
% y_DetShift=p(7).*cosd(tth/2)./tand(tth/2);
% 
% y=y_stress+y_VSShift+y_DetShift;
% end


function y=ahkl_ResultsEvaluation_advanced(pT,ss,GlobOpt,FitOpt,AdvancedShiftOptions)
y.a0=pT(1);
y.a0_err=ss(1);

if FitOpt.FitSF==1
    % p(2)=a0*alpha
    y.alphaSF_ef=pT(2);
    if FitOpt.Material.av==1
       y.alphaSF=pT(2)./y.a0;
    end
end

if GlobOpt.FitVerticalShift==1
    % p(3)=a0*sycos=a0*s/2R
    y.sycos_ef=pT(3);
    y.sycos_ef_err=ss(3);
    if FitOpt.Material.av==1
       y.sycos=pT(3)./y.a0;
       y.sycos_err=ss(3)./y.a0;
    end
end

if GlobOpt.FitDetectorShift==1
    % p(4)=a0*dth
    y.dth_ef=pT(4);
    y.dth_ef_err=pT(4);
    if FitOpt.Material.av==1
       y.dth=pT(4)./y.a0;
       y.dth_err=ss(4)./y.a0;
    end
end

FitN=AdvancedShiftOptions.FitN;
y.sigP=pT(5:4+FitN);
y.sigP_err=ss(5:4+FitN);

end


% function y=dhkl_ResultsEvaluation_hexagonal(pT,ss,x,ahkl,hkl,GlobOpt)
% h=hkl(:,1);
% k=hkl(:,2);
% l=hkl(:,3);
% switch GlobOpt.PeakShiftModel
%     case 'NoStress'
%         % p(1)=a0, p(2)=a0/c0
%         y.a0=pT(1);
%         y.a0_err=ss(1);
%         y.acrat=pT(2);
%         y.acrat_err=ss(2);
%         y.c0=pT(1)/pT(2);
%         y.c0_err=sqrt((ss(1)/pT(2))^2+(pT(1)*ss(2)/(pT(2))^2)^2);
%     case 'IsotropicStress'
%         y.a_norm=pT(4);
%         y.a_norm_err=ss(4);
%         y.a_par=pT(3)+pT(4);
%         y.a_par_err=ss(3)+ss(4);
%         % a0*(sig*(1/2*s2*sin2psi+2*s1)+1)
%         y.acrat=pT(2);
%         y.acrat_err=ss(2);
%         y.sig_ef=pT(3); % a0*sig*(s11+s12-2*s13)
%         y.sig_ef_err=ss(3); % a0*sig*(s11+s12-2*s13)
%         y.s13term=pT(4); % 2*a0*sig*s13+1
%         y.s13term_err=ss(4); % 2*a0*sig*s13+1
%         if GlobOpt.Material.av==1
%             y.sin2psi0=2*GlobOpt.Material.nu/(1+GlobOpt.Material.nu);
% 
%             y.a0=(y.a_par-y.a_norm)*y.sin2psi0+y.a_norm;
%             y.a0_err=(y.a_par_err-y.a_norm_err)*y.sin2psi0+y.a_norm_err;
%             y.sig=y.sig_ef*GlobOpt.Material.E/(y.a0.*(1+GlobOpt.Material.nu));
%             y.sig_err=y.sig_ef_err*GlobOpt.Material.E/(y.a0.*(1+GlobOpt.Material.nu))+y.sig_ef*y.a0_err*GlobOpt.Material.E/(y.a0.^2.*(1+GlobOpt.Material.nu));
%             y.c0=y.a0/y.acrat;
%             y.c0_err=sqrt((y.a0_err/y.acrat)^2+(y.a0*y.acrat_err/(y.acrat)^2)^2);
%         end
%     case 'VookWitt'
% %         % first the anisotropy needs to be substracted
% %         y.a_norm=pT(4);
% %         y.a_norm_err=ss(4);
% %         y.a_par=pT(1)+pT(4);
% %         y.a_par_err=ss(1)+ss(4);
% %         y.sig_ef=pT(1); % a0*sig*(1+nu_h00)/E_h00, but the stress needs to be evaluated from p(1) from the original fit, for cubic crystals E_h00=E and nu_h00=nu
% %         y.sig_ef_err=ss(1); % a0*sig*(1+nu_h00)/E_h00, but the stress needs to be evaluated from p(1) from the original fit, for cubic crystals E_h00=E and nu_h00=nu
% %         if GlobOpt.Material.av==1
% %             if isnan(GlobOpt.Material.E_100)
% %                 y.sin2psi0=2*GlobOpt.Material.nu/(1+GlobOpt.Material.nu);
% % 
% %                 y.a0=(y.a_par-y.a_norm)*y.sin2psi0+y.a_norm;
% %                 y.a0_err=(y.a_par_err-y.a_norm_err)*y.sin2psi0+y.a_norm_err;
% %                 y.sig=y.sig_ef*GlobOpt.Material.E/(y.a0.*(1+GlobOpt.Material.nu));
% %                 y.sig_err=y.sig_ef_err*GlobOpt.Material.E/(y.a0.*(1+GlobOpt.Material.nu))+y.sig_ef*y.a0_err*GlobOpt.Material.E/(y.a0.^2.*(1+GlobOpt.Material.nu));
% %             else
% %                 y.sin2psi0=2*GlobOpt.Material.nu_100/(1+GlobOpt.Material.nu_100);
% % 
% %                 y.a0=(y.a_par-y.a_norm)*y.sin2psi0+y.a_norm;
% %                 y.a0_err=(y.a_par_err-y.a_norm_err)*y.sin2psi0+y.a_norm_err;
% %                 y.sig=y.sig_ef*GlobOpt.Material.E_100/(y.a0.*(1+GlobOpt.Material.nu_100));
% %                 y.sig_err=y.sig_ef_err*GlobOpt.Material.E_100/(y.a0.*(1+GlobOpt.Material.nu_100))+y.a0_err*GlobOpt.Material.E_100/(y.a0.^2.*(1+GlobOpt.Material.nu_100));
% %             end
% %         end
%     case 'Reuss'
% %         % first the anisotropy needs to be substracted
% %         y.a_norm=pT(3);
% %         y.a_norm_err=ss(3);
% %         y.a_par=pT(1)+pT(3);
% %         y.a_par_err=ss(1)+ss(3);
% %         y.sig_ef=pT(1); % a0*sig*(1+nu_h00)/E_h00, but the stress needs to be evaluated from p(1) from the original fit, for cubic crystals E_h00=E and nu_h00=nu
% %         y.sig_ef_err=ss(1); % a0*sig*(1+nu_h00)/E_h00, but the stress needs to be evaluated from p(1) from the original fit, for cubic crystals E_h00=E and nu_h00=nu
% %         if GlobOpt.Material.av==1
% %             if isnan(GlobOpt.Material.E_100)
% %                 y.sin2psi0=2*GlobOpt.Material.nu/(1+GlobOpt.Material.nu);
% % 
% %                 y.a0=(y.a_par-y.a_norm)*y.sin2psi0+y.a_norm;
% %                 y.a0_err=(y.a_par_err-y.a_norm_err)*y.sin2psi0+y.a_norm_err;
% %                 y.sig=y.sig_ef*GlobOpt.Material.E/(y.a0.*(1+GlobOpt.Material.nu));
% %                 y.sig_err=y.sig_ef_err*GlobOpt.Material.E/(y.a0.*(1+GlobOpt.Material.nu))+y.sig_ef*y.a0_err*GlobOpt.Material.E/(y.a0.^2.*(1+GlobOpt.Material.nu));
% %             else
% %                 y.sin2psi0=2*GlobOpt.Material.nu_100/(1+GlobOpt.Material.nu_100);
% % 
% %                 y.a0=(y.a_par-y.a_norm)*y.sin2psi0+y.a_norm;
% %                 y.a0_err=(y.a_par_err-y.a_norm_err)*y.sin2psi0+y.a_norm_err;
% %                 y.sig=y.sig_ef*GlobOpt.Material.E_100/(y.a0.*(1+GlobOpt.Material.nu_100));
% %                 y.sig_err=y.sig_ef_err*GlobOpt.Material.E_100/(y.a0.*(1+GlobOpt.Material.nu_100))+y.a0_err*GlobOpt.Material.E_100/(y.a0.^2.*(1+GlobOpt.Material.nu_100));
% %             end
% %         end
% end
% 
% % if GlobOpt.FitSF==1
% %     % p(5)=a0*alpha
% %     y.alphaSF_ef=pT(5);
% %     if GlobOpt.Material.av==1
% %        y.alphaSF=pT(5)./y.a0;
% %     end
% % end
% 
% if GlobOpt.FitVerticalShift==1
%     % p(5)=a0*dth
%     y.sycos_ef=pT(5);
%     if GlobOpt.Material.av==1
%        y.sycos=pT(5)./y.a0;
%     end
% end
% 
% if GlobOpt.FitDetectorShift==1
%     % p(6)=a0*sycos=a0*s/2R
%     y.dth_ef=pT(6);
%     if GlobOpt.Material.av==1
%        y.dth=pT(6)./y.a0;
%     end
% end
% end


function [results,ahkl_T,iserr]=EvaluatePeakShiftAdvanced(F,ahkl,sin2psi,hkl,tth,psi,phi,GlobOpt,FitOpt,AdvancedShiftOptions)
% it is neccessary to set the parameters vector to fit
% 4 is for a0, stacking faults influence, vertical shift and detector shift
% corrections, the rest is the number of parameters for the stress fit
SigmaToFit_N=AdvancedShiftOptions.FitN;
p0=zeros(1,4+SigmaToFit_N);
fix=ones(1,4+SigmaToFit_N);
lb=-Inf*ones(1,4+SigmaToFit_N);
ub=+Inf*ones(1,4+SigmaToFit_N);

p0(1)=mean(ahkl);
fix(1)=0;
lb(1)=0.8*p0(1);
ub(1)=1.2*p0(1);

if FitOpt.FitSF==1
    if strcmp(FitOpt.sys,'cubic_fcc') % for bcc and hcp there is no influence on peak shift
        p0(2)=1e-2; % it is zero by default
        fix(2)=0;
        lb(2)=0;
        ub(2)=1;
    end
end

if GlobOpt.FitVerticalShift==1
    p0(3)=0;
    fix(3)=0;
    lb(3)=-0.1; % in radians
    ub(3)=0.1; % in radians
end

if GlobOpt.FitDetectorShift==1
    p0(4)=0;
    fix(4)=0;
    lb(4)=-0.1; % in radians
    ub(4)=0.1; % in radians
end

pT0=polyfit(sin2psi,ahkl,1);
sighlp=pT0(1)/mean(F,'all')/p0(1);
% p0(5:4+SigmaToFit_N)=sighlp;
fix(5:4+SigmaToFit_N)=0;
lb(5:4+SigmaToFit_N)=-50*abs(sighlp);
ub(5:4+SigmaToFit_N)=50*abs(sighlp);

if numel(tth)<=sum(fix==0)
    iserr=1; % more parameters than datapoints
elseif sum(imag(p0))~=0
    iserr=2; % initial guess gives complex parameters
else
    iserr=0;
    [pT1,resnorm,~,~,~,~,jacob]=lsqcurvefit(@(p,xdata) ahkl_shift_advanced(xdata,interlace(p0,p,fix),F,hkl,tth,psi,phi,GlobOpt,FitOpt,AdvancedShiftOptions),p0(~fix),sin2psi,ahkl,lb(~fix),ub(~fix));
    ss1=full(diag(sqrt(resnorm*inv(jacob'*jacob))/(length(sin2psi)-length(pT1))));
    pT=interlace(p0,pT1,fix); % for the correct evaluation it is necessary to have all parameters, included those fixed to zero  
    ss=interlace(zeros(size(pT)),ss1,fix);

    if sum(imag([pT1(:); ss1(:)]))~=0 || max(abs(ss1(:)./pT1(:)))>0.5
        iserr=-1; % complex values of fitting parameters, or large errors
    end
end

if iserr<=0
    results=ahkl_ResultsEvaluation_advanced(pT,ss,GlobOpt,FitOpt,AdvancedShiftOptions);
    ahkl_T=ahkl_shift_advanced(sin2psi,pT,F,hkl,tth,psi,phi,GlobOpt,FitOpt,AdvancedShiftOptions);
else
    results=NaN;
    ahkl_T=NaN;
end

end


% function [results,dhkl_T,iserr]=EvaluatePeakShift_hexagonal(sin2psi,dhkl,tth,hkl,GlobOpt)
% p0=zeros(1,7);
% fix=ones(1,7);
% lb=-Inf*ones(1,7);
% ub=+Inf*ones(1,7);
% [a0,c0]=EstimateHexagonalLattice(hkl,dhkl);
% switch GlobOpt.PeakShiftModel
%     case 'NoStress'
%         p0(1:2)=[a0 a0/c0];
%         fix(1:4)=[0 0 1 1];
%         lb(1:2)=[0.8*p0(1) 0];
%         ub(1:2)=[1.2*p0(1) p0(2)*10];
%     case 'IsotropicStress'
%         p0(1:2)=[a0 a0/c0];
%         lb(1:2)=[0.8*p0(1) 0];
%         ub(1:2)=[1.2*p0(1) p0(2)*10];
%         if GlobOpt.FitHcpLP==1
%             fix(1:4)=[1 0 0 0];
%         else
%             fix(1:4)=[1 1 0 0];
%         end
%     case 'VookWitt'
%         disp('Sorry, I can''t handle this method yet.')
% %         pT0=polyfit(sin2psi,ahkl,1);
% %         p0(1:2)=pT0;
% %         fix(1:4)=[0 0 0 0];
%     case 'Reuss'
%         disp('Sorry, I can''t handle this method yet.')
% %         pT0=polyfit(sin2psi,ahkl,1);
% %         p0(1:2)=pT0;
% %         fix(1:4)=[0 0 0 1];
% end
% 
% % % There is no peak shift
% % if GlobOpt.FitSF==1
% %     p0(6)=1e-2; % it is zero by default
% %     fix(6)=0;
% %     lb(6)=0;
% %     ub(6)=1; % is the alpha density or probability?
% % end
% 
% if GlobOpt.FitVerticalShift==1
%     p0(5)=0;
%     fix(5)=0;
%     lb(5)=-0.1; % in radians
%     ub(5)=0.1; % in radians
% end
% 
% if GlobOpt.FitDetectorShift==1
%     p0(6)=0;
%     fix(6)=0;
%     lb(6)=-0.1; % in radians
%     ub(6)=0.1; % in radians
% end
% 
% if numel(tth)<=sum(fix==0)
%     iserr=1; % more parameters than datapoints
% elseif sum(imag(p0))~=0
%     iserr=2; % initial guess gives complex parameters
% else
%     iserr=0;
%     [pT1,resnorm,~,~,~,~,jacob]=lsqcurvefit(@(p,xdata) dhkl_shift_hexagonal(xdata,interlace(p0,p,fix),hkl,tth,GlobOpt),p0(~fix),sin2psi,dhkl,lb(~fix),ub(~fix));
%     ss1=full(diag(sqrt(resnorm*inv(jacob'*jacob))/(length(sin2psi)-length(pT1))));
%     pT=interlace(p0,pT1,fix); % for the correct evaluation it is necessary to have all parameters, included those fixed to zero
%     ss=interlace(zeros(size(pT)),ss1,fix);
% 
%     if sum(imag([pT1(:); ss1(:)]))~=0 || max(abs(ss1(:)./pT1(:)))>0.5
%         iserr=-1; % complex values of fitting parameters, or large errors
%     end
% end
% 
% if iserr<=0
%     results=dhkl_ResultsEvaluation_hexagonal(pT,ss,sin2psi,dhkl,hkl,GlobOpt);
%     dhkl_T=dhkl_shift_hexagonal(sin2psi,pT,hkl,tth,GlobOpt);
% else
%     results=NaN;
%     dhkl_T=NaN;
% end
% 
% end