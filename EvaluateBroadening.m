function [results,broad_0,broad_T,iserr]=EvaluateBroadening(broaddata,data,GlobOpt,FitOpt)
% The reason of input broaddata separately is that this function is called
% separately to evaluate broadening from FWHM and integral width
% (IntegratedIntensity/Imax)
hkl=[data.h data.k data.l];
tth=data.tth;
psi=data.psi;
phi=data.phi;
% The crucial part here is to correctly recompute the angular broadening
% into reciprocal space.
broad_0=1/GlobOpt.lam*cosd(tth/2).*broaddata*pi/360; % I need to recompute it to the rms of q-vector, fwhm is measured in 2th, therefore also devide it by 2, there is no 4*pi to be in correspondence with Scherrer equation
[results,broad_T,iserr]=EvaluatePeakBroadening(tth,psi,phi,broad_0,hkl,GlobOpt,FitOpt);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% AUXILIARY FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=broad_hkl(tth,p,psi,phi,hkl,GlobOpt,FitOpt)
% Kernel function for the broadening fit.
% Depending on the model for the peaks broadening, different number
% of fitting parameters are needed (p(1:4)). Therefore, some of them need
% to be fixed to 0 during the fit depending on the choosen model.
% The evaluation of dislocation density, SF probability and other constants
% is done separately after fit, in dependence on the model and also
% in dependence on availability of material constants input by the user.
h=hkl(:,1);
k=hkl(:,2);
l=hkl(:,3);

% Size
% p(1)=0.94/D
y_size=zeros(size(tth));
if FitOpt.FitSize==1
    switch FitOpt.SizeOptions.GrownDirectionOpt
        case 'Sample normal'
            psi_hkl=psi;
            phi_hkl=phi-FitOpt.SizeOptions.GrownDirectionPerpAzimuth;
        case 'Specific crystallographic direction'
            DomDir=FitOpt.SizeOptions.GrownDirectionHKL;
            DomDir2=FitOpt.SizeOptions.GrownDirectionHKL2;
            psi_hkl=acosd((hkl(:,1)*DomDir(1)+hkl(:,2)*DomDir(2)+hkl(:,3)*DomDir(3))./(vecnorm(hkl')'.*norm(DomDir)));
            % Projection of hkl onto plane perpendicular to HKL
            hkl_par=hkl-(hkl(:,1)*DomDir(1)+hkl(:,2)*DomDir(2)+hkl(:,3)*DomDir(3))/dot(DomDir,DomDir).*DomDir;
            phi_hkl=acosd((hkl_par(:,1)*DomDir2(1)+hkl_par(:,2)*DomDir2(2)+hkl_par(:,3)*DomDir2(3))./(vecnorm(hkl_par')'.*norm(DomDir)));
            phi_hkl(isnan(phi_hkl))=0;
        otherwise
            disp('Uknown option for the grown direction.');
            psi_hkl=NaN;
            phi_hkl=NaN;
    end

    switch FitOpt.GrainsShape
        case 'Sphere'
            y_size=0.94/p(1);
        case 'Cylinder'
            D=p(1)./cosd(psi_hkl);
            psi0=atand(p(2)/p(1));
            ind=find(psi_hkl>psi0);
            D(ind)=p(2)./sind(psi_hkl(ind));
            y_size=0.94./D;
        case 'Cylinder with elliptical base'
            D=p(1)./cosd(psi_hkl);
            Dpar=p(2).*p(3)./sqrt(p(2).^2.*sind(phi_hkl).^2+p(3).^2.*cosd(phi_hkl).^2);
            psi0=atand(Dpar./p(1));
            ind=find(psi_hkl>psi0);
            D(ind)=Dpar(ind)./sind(psi_hkl(ind));
            y_size=0.94./D;
        case 'Rotational ellipsoid'
            D=p(1)*p(2)./sqrt(p(1).^2.*sind(psi_hkl).^2+p(2).^2.*cosd(psi_hkl).^2);
            y_size=0.94./D;
        case 'General ellipsoid'
            D=p(1)*p(2)*p(3)./sqrt(p(1).^2*(p(3).^2.*cosd(phi_hkl).^2+p(2).^2.*sind(phi_hkl).^2).*sind(psi_hkl).^2+p(2).^2*p(3).^2.*cosd(psi_hkl).^2);
            y_size=0.94./D;
    end
end

% Stacking faults
y_SF=zeros(size(tth));
if FitOpt.FitSF==1
    switch FitOpt.sys
        case 'cubic_fcc'
            % p(2)/a*sum_affected(abs(h+k+l)/(m_hkl*sqrt(h^2+k^2+l^2)));
            for n=1:numel(h)
                [planes,m_hkl]=GetEquivalentPlanes([h(n) k(n) l(n)],'cubic');
                L=planes(:,1)+planes(:,2)+planes(:,3); % h+k+l
                % lines with mod(L,3)==0 are not affected
                b=ones(size(L));
                b(mod(L,3)==0)=0;
                y_SF(n)=y_SF(n)+p(4).*sum(b.*abs(planes(:,1)+planes(:,2)+planes(:,3))./(m_hkl.*sqrt(planes(:,1).^2+planes(:,2).^2+planes(:,3).^2)));
            end
        case 'cubic_bcc'
            % p(2)/a*sum_affected(abs(-h-k+2*l)/(m_hkl*sqrt(h^2+k^2+l^2)));
            for n=1:numel(h)
                [planes,m_hkl]=GetEquivalentPlanes([h(n) k(n) l(n)],'cubic');
                L=-planes(:,1)-planes(:,2)+2*planes(:,3); % -h-k+2*l
                % lines with mod(L,3)==0 are not affected
                b=ones(size(L));
                b(mod(L,3)==0)=0;
%                 b(mod(L,3)==2)=-1;
                y_SF(n)=y_SF(n)+p(4).*sum(b.*abs(-planes(:,1)-planes(:,2)+2*planes(:,3))./(m_hkl.*sqrt(planes(:,1).^2+planes(:,2).^2+planes(:,3).^2)));
            end
        case 'hexagonal'
            % 3*p(2)*l*dhkl/c^2
            dhkl=GlobOpt.lam./(2*sind(tth/2));
            [a0,c0]=EstimateHexagonalLattice(hkl,dhkl);
            B=Reci(a0,a0,c0,90,90,120);
            dhkl0=zeros(size(tth));
            for n=1:numel(tth)
               dhkl0(n)=2*pi./norm(B*hkl(n,:)');
            end
            y_SF=3*p(4).*hkl(:,3).*dhkl0/c0^2;
    end
end

% Microstrain
switch FitOpt.MicrostrainOptions.MicrostrainModel
    case 'none'
        y_strain=zeros(size(tth));
    case 'Isotropic'
        % p(3)*sin(th)
        y_strain=p(5).*sind(tth/2);
    case 'Dislocations (polycrystall)'
        switch FitOpt.sys
            case {'cubic_fcc'; 'cubic_bcc'}
                % p(2)*sqrt(1-p(3)*Gamma)*sin(th)
                Gamma=(h.^2.*k.^2+k.^2.*l.^2+h.^2.*l.^2)./((h.^2+k.^2+l.^2).^2);
                y_strain=p(5).*sqrt(1-p(6).*Gamma).*sind(tth/2);
            case 'hexagonal'
                % p(2)*sqrt(p(3)*ds_hkl.^2+((2*p(4)*(hkl(:,1).^2+hkl(:,1).*hkl(:,2)+hkl(:,2).^2)+p(5).*hkl(:,3).^2).*hkl(:,3).^2)./(9/4*(a0.*ds_hkl).^2))
                dhkl=GlobOpt.lam./(2*sind(tth/2));
                [a0,c0]=EstimateHexagonalLattice(hkl,dhkl);
                B=Reci(a0,a0,c0,90,90,120);
                dhkl0=zeros(size(tth));
                for n=1:numel(tth)
                    dhkl0(n)=2*pi./norm(B*hkl(n,:)');
                end
                ds_hkl=1./dhkl0;
                y_strain=p(5)*sqrt(p(6)*ds_hkl.^2+((2*p(7)*(hkl(:,1).^2+hkl(:,1).*hkl(:,2)+hkl(:,2).^2)+p(8).*hkl(:,3).^2).*hkl(:,3).^2)./(9/4*(a0.*ds_hkl).^2));
        end
    case 'Dislocations (specific Burgers vector)'
        eq=[sind(psi).*cosd(phi) sind(psi).*sind(phi) cosd(psi)]; % list of vectors
        eB=[sind(p(6)).*cosd(p(7)) sind(p(6)).*sind(p(7)) cosd(p(6))]; % one vector
        sp=eq(:,1)*eB(1)+eq(:,2)*eB(2)+eq(:,3)*eB(:,3);
        y_strain=p(5).*abs(sp);
end

% in the end, size and strain contributions are sum with the NFitOrder
N=FitOpt.BroadeningNFitOrder;
y=((y_size+y_SF).^N+y_strain.^N).^(1/N);
end


function y=broadening_ResultsEvaluation(pT,ss,GlobOpt,FitOpt)
% Size
if FitOpt.FitSize==1
    switch FitOpt.GrainsShape
        case 'Sphere'
            y.D=pT(1);
            y.D_err=ss(1);
        case 'Cylinder'
            y.D_norm=pT(1);
            y.D_norm_err=ss(1);
            y.D_par=pT(2);
            y.D_par_err=ss(2);
        case 'Cylinder with elliptical base'
            y.D1=pT(1);
            y.D1_err=ss(1);
            y.D2=pT(2);
            y.D2_err=ss(2);
            y.D3=pT(3);
            y.D3_err=ss(3);
        case 'Rotational ellipsoid'
            y.D_norm=pT(1);
            y.D_norm_err=ss(1);
            y.D_par=pT(2);
            y.D_par_err=ss(2);
        case 'General ellipsoid'
            y.D1=pT(1);
            y.D1_err=ss(1);
            y.D2=pT(2);
            y.D2_err=ss(2);
            y.D3=pT(3);
            y.D3_err=ss(3);
    end
end

% Stacking faults
if FitOpt.FitSF==1
    y.al_ef=pT(4);
    y.al_ef_err=ss(4);
end

% Microstrain
% switch FitOpt.MicrostrainOptions.PeakBroadeningModel
switch FitOpt.MicrostrainOptions.MicrostrainModel
    case 'none'
        % no microstrain parameters fitted
    case 'Isotropic'
        % broad_strain=4*e/lam*sin(th)
        y.e=pT(5)*GlobOpt.lam/4;
        y.e_err=ss(5)*GlobOpt.lam/4;
    case 'Dislocations (polycrystall)'
        switch FitOpt.sys
            case {'cubic_fcc'; 'cubic_bcc'} 
                % broad_strain_dis=4/lam*e_100*sqrt(1-q*Gamma)*sin(th)      
                y.e_dis_100=pT(5)*GlobOpt.lam/4;
                y.e_dis_100_err=ss(5)*GlobOpt.lam/4;
                y.q_dis=pT(6);
                y.q_dis_err=ss(6);
            case 'hexagonal'
                % p(2)*sqrt(p(3)*ds_hkl.^2+((2*p(4)*(hkl(:,1).^2+hkl(:,1).*hkl(:,2)+hkl(:,2).^2)+p(5).*hkl(:,3).^2).*hkl(:,3).^2)./(9/4*(a0.*ds_hkl).^2))
                y.e_dis_100=pT(5)*GlobOpt.lam/4;
                y.e_dis_100_err=ss(5)*GlobOpt.lam/4;
                y.al_dis=pT(6);
                y.al_dis_err=ss(6);
                y.be_dis=pT(7);
                y.be_dis_err=ss(7);
                y.ga_dis=pT(8);
                y.ga_dis_err=ss(8);
        end
    case 'Dislocations (specific Burgers vector)'
         y.e_dis_100=pT(5)*GlobOpt.lam/4;
         y.e_dis_100_err=ss(5)*GlobOpt.lam/4;
         y.PsiB=pT(6);
         y.PsiB_err=ss(6);
         y.PhiB=pT(7);
         y.PhiB_err=ss(7);
end
end

function [results,broad_T,iserr]=EvaluatePeakBroadening(tth,psi,phi,broad,hkl,GlobOpt,FitOpt)
% p0=FitOpt.SizeOptions.ParGuess,FitOpt.SFOptions.ParGuess,FitOpt.MicrostrainOptions.ParGuess
p0=zeros(1,8); % D, al_SF_ef, isostrain/k1_dis, k2_dis, k3_dis, k4_dis
fix=ones(1,8);
lb=-Inf*ones(1,8);
ub=+Inf*ones(1,8);
N=FitOpt.BroadeningNFitOrder;
% switch FitOpt.MicrostrainOptions.PeakBroadeningModel
switch FitOpt.MicrostrainOptions.MicrostrainModel
    % 'none', 'Isotropic', 'Dislcoations'
    case 'none'
        % p(1)=0.94/D;
        p0(1)=0.94/mean(broad);
        fix(1)=0;
        lb(1)=p0(1)/100;
        ub(1)=p0(1)*100;
    case 'Isotropic'
        % 0.94/D+4*e/lambda*sin(th)
        pT0=polyfit((sind(tth/2)).^N,broad.^N,1);
        p0(1)=0.94/max([pT0(2) mean(broad.^N)]).^(1/N);
        p0(5)=max([pT0(1).^(1/N)*GlobOpt.lam/4 0.005]);
        fix([1 5])=[0 0];
        lb([1 5])=[p0(1)/100 p0(5)/100];
        ub([1 5])=[p0(1)*100 p0(5)*100];
    case 'Dislocations (polycrystall)'
        switch FitOpt.sys
            case {'cubic_fcc'; 'cubic_bcc'} 
                fix([1 5 6])=[0 0 0];
                pT0=polyfit((sind(tth/2)).^N,broad.^N,1);
                p0(1)=0.94/max([pT0(2) mean(broad.^N)]).^(1/N);
                lb(1)=p0(1)/100;
                ub(1)=p0(1)*100;

                % p(2)*sqrt(1-p(3)*Gamma)*sin(th)
                p0(5)=max([pT0(1).^(1/N) 0.001]);
                p0(6)=0;
                lb(5:6)=[p0(5)/100 -100];
                ub(5:6)=[p0(5)*100 100];
            case 'hexagonal'
                fix([1 5 6 7])=[0 0 0 0 0];
                pT0=polyfit((sind(tth/2)).^N,broad.^N,1);
                p0(1)=0.94/max([pT0(2) mean(broad.^N)]).^(1/N);
                lb(1)=p0(1)/100;
                ub(1)=p0(1)*100;

                % p(2)*sqrt(1-p(3)*Gamma)*sin(th)
                p0(5)=max([pT0(1).^(1/N) 0.001]);
                p0(6:8)=[0 0 0];
                lb(5:8)=[p0(5)/100 -100 -100 -100];
                ub(5:8)=[p0(5)*100 100 100 100];
        end
    case 'Dislocations (specific Burgers vector)'
          fix([5 6 7])=[0 0 0];
          pT0=polyfit((sind(tth/2)).^N,broad.^N,1);
          p0(1)=0.94/max([pT0(2) mean(broad.^N)]).^(1/N);
          lb(1)=p0(1)/100;
          ub(1)=p0(1)*100;
          % p(2)*sqrt(1-p(3)*Gamma)*sin(th)
          p0(5)=max([pT0(1).^(1/N) 0.001]);
          p0(6:7)=[0 0];
          lb(5:7)=[p0(5)/100 -90 -90];
          ub(5:7)=[p0(5)*100 90 90];
end

% Grain size and shape
p0(2)=p0(1); p0(3)=p0(1);
lb(2)=lb(1); lb(3)=lb(1);
ub(2)=ub(1); ub(3)=ub(1);
if FitOpt.FitSize
    switch FitOpt.GrainsShape
        case 'Sphere'
            fix(1:3)=[0 1 1];
            % if there is a guess value, use it
            if ~isnan(FitOpt.SizeOptions.ParGuess)
                p0(1)=FitOpt.SizeOptions.ParGuess;
                if FitOpt.SizeOptions.ParFix
                    fix(1)=1;
                end
            end
            if ~isnan(FitOpt.SizeOptions.ParLB)
                lb(1)=FitOpt.SizeOptions.ParLB;
            end
            if ~isnan(FitOpt.SizeOptions.ParUB)
                ub(1)=FitOpt.SizeOptions.ParUB;
            end
        case 'Cylinder'
            fix(1:3)=[0 0 1];
            % if there is a guess value, use it
            if ~isnan(FitOpt.SizeOptions.ParGuess(1))
                p0(1)=FitOpt.SizeOptions.ParGuess(1);
                if FitOpt.SizeOptions.ParFix(1)
                    fix(1)=1;
                end
            end
            if ~isnan(FitOpt.SizeOptions.ParLB(1))
                lb(1)=FitOpt.SizeOptions.ParLB(1);
            end
            if ~isnan(FitOpt.SizeOptions.ParUB(1))
                ub(1)=FitOpt.SizeOptions.ParUB(1);
            end
            % if there is a guess value, use it
            if ~isnan(FitOpt.SizeOptions.ParGuess(2))
                p0(2)=FitOpt.SizeOptions.ParGuess(2);
                if FitOpt.SizeOptions.ParFix(2)
                    fix(2)=1;
                end
            end
            if ~isnan(FitOpt.SizeOptions.ParLB(2))
                lb(2)=FitOpt.SizeOptions.ParLB(2);
            end
            if ~isnan(FitOpt.SizeOptions.ParUB(2))
                ub(2)=FitOpt.SizeOptions.ParUB(2);
            end
        case 'Cylinder with elliptical base'
            fix(1:3)=[0 0 0];
            % if there is a guess value, use it
            if ~isnan(FitOpt.SizeOptions.ParGuess(1))
                p0(1)=FitOpt.SizeOptions.ParGuess(1);
                if FitOpt.SizeOptions.ParFix(1)
                    fix(1)=1;
                end
            end
            if ~isnan(FitOpt.SizeOptions.ParLB(1))
                lb(1)=FitOpt.SizeOptions.ParLB(1);
            end
            if ~isnan(FitOpt.SizeOptions.ParUB(1))
                ub(1)=FitOpt.SizeOptions.ParUB(1);
            end
            % if there is a guess value, use it
            if ~isnan(FitOpt.SizeOptions.ParGuess(2))
                p0(2)=FitOpt.SizeOptions.ParGuess(2);
                if FitOpt.SizeOptions.ParFix(2)
                    fix(2)=1;
                end
            end
            if ~isnan(FitOpt.SizeOptions.ParLB(2))
                lb(2)=FitOpt.SizeOptions.ParLB(2);
            end
            if ~isnan(FitOpt.SizeOptions.ParUB(2))
                ub(2)=FitOpt.SizeOptions.ParUB(2);
            end
            % if there is a guess value, use it
            if ~isnan(FitOpt.SizeOptions.ParGuess(3))
                p0(3)=FitOpt.SizeOptions.ParGuess(3);
                if FitOpt.SizeOptions.ParFix(3)
                    fix(3)=1;
                end
            end
            if ~isnan(FitOpt.SizeOptions.ParLB(3))
                lb(3)=FitOpt.SizeOptions.ParLB(3);
            end
            if ~isnan(FitOpt.SizeOptions.ParUB(3))
                ub(3)=FitOpt.SizeOptions.ParUB(3);
            end
        case 'Rotational ellipsoid'
            fix(1:3)=[0 0 1];
            % if there is a guess value, use it
            if ~isnan(FitOpt.SizeOptions.ParGuess(1))
                p0(1)=FitOpt.SizeOptions.ParGuess(1);
                if FitOpt.SizeOptions.ParFix(1)
                    fix(1)=1;
                end
            end
            if ~isnan(FitOpt.SizeOptions.ParLB(1))
                lb(1)=FitOpt.SizeOptions.ParLB(1);
            end
            if ~isnan(FitOpt.SizeOptions.ParUB(1))
                ub(1)=FitOpt.SizeOptions.ParUB(1);
            end
            % if there is a guess value, use it
            if ~isnan(FitOpt.SizeOptions.ParGuess(2))
                p0(2)=FitOpt.SizeOptions.ParGuess(2);
                if FitOpt.SizeOptions.ParFix(2)
                    fix(2)=1;
                end
            end
            if ~isnan(FitOpt.SizeOptions.ParLB(2))
                lb(2)=FitOpt.SizeOptions.ParLB(2);
            end
            if ~isnan(FitOpt.SizeOptions.ParUB(2))
                ub(2)=FitOpt.SizeOptions.ParUB(2);
            end
        case 'General ellipsoid'
            fix(1:3)=[0 0 0];
            % if there is a guess value, use it
            if ~isnan(FitOpt.SizeOptions.ParGuess(1))
                p0(1)=FitOpt.SizeOptions.ParGuess(1);
                if FitOpt.SizeOptions.ParFix(1)
                    fix(1)=1;
                end
            end
            if ~isnan(FitOpt.SizeOptions.ParLB(1))
                lb(1)=FitOpt.SizeOptions.ParLB(1);
            end
            if ~isnan(FitOpt.SizeOptions.ParUB(1))
                ub(1)=FitOpt.SizeOptions.ParUB(1);
            end
            % if there is a guess value, use it
            if ~isnan(FitOpt.SizeOptions.ParGuess(2))
                p0(2)=FitOpt.SizeOptions.ParGuess(2);
                if FitOpt.SizeOptions.ParFix(2)
                    fix(2)=1;
                end
            end
            if ~isnan(FitOpt.SizeOptions.ParLB(2))
                lb(2)=FitOpt.SizeOptions.ParLB(2);
            end
            if ~isnan(FitOpt.SizeOptions.ParUB(2))
                ub(2)=FitOpt.SizeOptions.ParUB(2);
            end
            % if there is a guess value, use it
            if ~isnan(FitOpt.SizeOptions.ParGuess(3))
                p0(3)=FitOpt.SizeOptions.ParGuess(3);
                if FitOpt.SizeOptions.ParFix(3)
                    fix(3)=1;
                end
            end
            if ~isnan(FitOpt.SizeOptions.ParLB(3))
                lb(3)=FitOpt.SizeOptions.ParLB(3);
            end
            if ~isnan(FitOpt.SizeOptions.ParUB(3))
                ub(3)=FitOpt.SizeOptions.ParUB(3);
            end
    end
else
    % if the sizes are not fitted, the grains are very big
    p0(1:2)=1e6;
    lb(1:2)=1e5;
    ub(1:2)=1e7;
    fix(1:2)=[1 1];
end

if FitOpt.FitSF==1
    % p(2)/a*sum_affected(abs(h+k+l)/(m_hkl*sqrt(h^2+k^2+l^2))); 'cubic_fcc'
    % p(2)/a*sum_affected(abs(-h-k+2*l)/(m_hkl*sqrt(h^2+k^2+l^2))); 'cubic_bcc'
    % 3*p(2)*l*dhkl/c^2; 'hexagonal'
    p0(4)=5e-3;
    fix(4)=0;
    lb(4)=0;
    ub(4)=1; % this is probability divided by lattice parameter, so the upper limit might be smaller
end

if FitOpt.FitMicrostrain
    switch FitOpt.MicrostrainOptions.MicrostrainModel
        case 'Isotropic'
            fix([5 6 7 8])=[0 1 1 1];
            % if there is a guess value, use it
            if ~isnan(FitOpt.MicrostrainOptions.ParGuess)
                p0(5)=FitOpt.MicrostrainOptions.ParGuess*4/GlobOpt.lam;
                if FitOpt.MicrostrainOptions.ParFix
                    fix(5)=1;
                end
            end
            if ~isnan(FitOpt.MicrostrainOptions.ParLB)
                lb(5)=FitOpt.MicrostrainOptions.ParLB*4/GlobOpt.lam;
            end
            if ~isnan(FitOpt.MicrostrainOptions.ParUB)
                ub(5)=FitOpt.MicrostrainOptions.ParUB*4/GlobOpt.lam;
            end
        case 'Dislocations (polycrystall)'
            fix([5 6 7 8])=[0 0 1 1];
            % if there is a guess value, use it
            if ~isnan(FitOpt.MicrostrainOptions.ParGuess(1))
                p0(5)=FitOpt.MicrostrainOptions.ParGuess(1)*4/GlobOpt.lam;
                if FitOpt.MicrostrainOptions.ParFix(1)
                    fix(5)=1;
                end
            end
            if ~isnan(FitOpt.MicrostrainOptions.ParLB(1))
                lb(5)=FitOpt.MicrostrainOptions.ParLB(1)*4/GlobOpt.lam;
            end
            if ~isnan(FitOpt.MicrostrainOptions.ParUB(1))
                ub(5)=FitOpt.MicrostrainOptions.ParUB(1)*4/GlobOpt.lam;
            end
            % if there is a guess value, use it
            if ~isnan(FitOpt.MicrostrainOptions.ParGuess(2))
                p0(6)=FitOpt.MicrostrainOptions.ParGuess(2);
                if FitOpt.MicrostrainOptions.ParFix(2)
                    fix(6)=1;
                end
            end
            if ~isnan(FitOpt.MicrostrainOptions.ParLB(2))
                lb(6)=FitOpt.MicrostrainOptions.ParLB(2);
            end
            if ~isnan(FitOpt.MicrostrainOptions.ParUB(2))
                ub(6)=FitOpt.MicrostrainOptions.ParUB(2);
            end
        case 'Dislocations (specific Burgers vector)'
            fix([5 6 7 8])=[0 0 0 1];
            % if there is a guess value, use it
            if ~isnan(FitOpt.MicrostrainOptions.ParGuess(1))
                p0(5)=FitOpt.MicrostrainOptions.ParGuess(1)*4/GlobOpt.lam;
                if FitOpt.MicrostrainOptions.ParFix(1)
                    fix(5)=1;
                end
            end
            if ~isnan(FitOpt.MicrostrainOptions.ParLB(1))
                lb(5)=FitOpt.MicrostrainOptions.ParLB(1)*4/GlobOpt.lam;
            end
            if ~isnan(FitOpt.MicrostrainOptions.ParUB(1))
                ub(5)=FitOpt.MicrostrainOptions.ParUB(1)*4/GlobOpt.lam;
            end
            % if there is a guess value, use it
            if ~isnan(FitOpt.MicrostrainOptions.ParGuess(2))
                p0(6)=FitOpt.MicrostrainOptions.ParGuess(2);
                if FitOpt.MicrostrainOptions.ParFix(2)
                    fix(6)=1;
                end
            end
            if ~isnan(FitOpt.MicrostrainOptions.ParLB(2))
                lb(6)=FitOpt.MicrostrainOptions.ParLB(2);
            end
            if ~isnan(FitOpt.MicrostrainOptions.ParUB(2))
                ub(6)=FitOpt.MicrostrainOptions.ParUB(2);
            end
            % if there is a guess value, use it
            if ~isnan(FitOpt.MicrostrainOptions.ParGuess(3))
                p0(7)=FitOpt.MicrostrainOptions.ParGuess(3);
                if FitOpt.MicrostrainOptions.ParFix(3)
                    fix(7)=1;
                end
            end
            if ~isnan(FitOpt.MicrostrainOptions.ParLB(3))
                lb(7)=FitOpt.MicrostrainOptions.ParLB(3);
            end
            if ~isnan(FitOpt.MicrostrainOptions.ParUB(3))
                ub(7)=FitOpt.MicrostrainOptions.ParUB(3);
            end
    end
end
if numel(tth)<=sum(fix==0)
    iserr=1; % more parameters than datapoints
elseif sum(imag(p0))~=0
    iserr=2; % initial guess gives complex parameters
else
    iserr=0;

    [pT1,resnorm,~,~,~,~,jacob]=lsqcurvefit(@(p,xdata) broad_hkl(xdata,interlace(p0,p,fix),psi,phi,hkl,GlobOpt,FitOpt),p0(~fix),tth,broad,lb(~fix),ub(~fix));
    ss1=full(diag(sqrt(resnorm*inv(jacob'*jacob))/(length(tth)-length(pT1))));
    pT=interlace(p0,pT1,fix); % for the correct evaluation it is necessary to have all parameters, included those fixed to zero
    ss=interlace(zeros(size(pT)),ss1,fix);

    if sum(imag([pT1(:); ss1(:)]))~=0 || max(abs(ss1(:)./pT1(:)))>0.5
        iserr=-1; % complex values of fitting parameters, or large errors
    end
end

if iserr<=0
    results=broadening_ResultsEvaluation(pT,ss,GlobOpt,FitOpt);
    broad_T=broad_hkl(tth,pT,psi,phi,hkl,GlobOpt,FitOpt);
else
    results=NaN;
    broad_T=NaN;
end

end