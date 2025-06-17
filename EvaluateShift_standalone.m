%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
label='Cr on MgO(110)';
filename='Fit results_Cr_on_MgO_110.dat';
GlobOpt.lam=1.5406; % CuKal1 wavelength
minint=1e2; % take only results with the intensity higher than minimal value
nu=0.21; % Poisson ratio
E=279; % Young modulus (GPa)
% Contrast factor obtained from anizc (http://metal.elte.hu/anizc)
% C100=0.178108196; q=1.126327936; b=2.497617265; % screw dislocation <111>
% C100=0.133443048; q=-0.119868455; b=2.039295957; % screw dislocation <110>
% C100=0.129114063; q=-2.001409881; b=2.497617265; % edge dislocation <111>{110}
C100=0.131904962; q=-1.849758323; b=2.497617265; % edge dislocation <111>{211}

% profile function parameters (Caglioti polynom)
u0=-0.0069; u1=0.0002; % U(chi)=u0+u1*chi
v0=0.0225; v1=-0.0002; % V(chi)=v0+v1*chi
w0=0.0572; w1=0.0001; % V(chi)=v0+v1*chi
Eta0=0.0624; Eta1=0.0719; % Eta=Eta0+Eta1*2Theta
Asym0=0.9924; Asym1=0.0416; % Asym=Asym0+Asym1/sin(2Theta)

GlobOpt.PeakShiftModel='Reuss'; % 'NoStress', 'Isotropic', 'VookWitt', 'Reuss'
GlobOpt.FitSF=1;
GlobOpt.FitShift=0;
% PeakBroadeningModel='betahkl_disl_cubic'; % betahkl_classic, betahkl_disl_cubic
GlobOpt.Material.av=1;
GlobOpt.Material.nu=0.21; % Poisson ratio
GlobOpt.Material.E=279; % Young modulus (GPa)
GlobOpt.sys='cubic_bcc'; % 'cubic_fcc', 'cubic_bcc', 'hexagonal'

%%%% DATA IMPORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts = detectImportOptions(filename);
opts.VariableNames={'h','k','l','grain','chi','phi','IntInt','tth','eta','FWHM_left','FWHM_right'}; % MatLab will read the hkl as a double, this will override it, however, one could overwrite hkl into h, k, l directly in the table
mytable=readtable(filename,opts);

psi=mytable.chi;
tth=mytable.tth;
hkl=[mytable.h(:) mytable.k(:) mytable.l(:)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PEAK SHIFT EVALUATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sin2psi=(sind(psi)).^2;
switch GlobOpt.sys
    case {'cubic_fcc','cubic_bcc'}
        ahkl=GlobOpt.lam./(2*sind(mytable.tth/2)).*sqrt(mytable.h.^2+mytable.k.^2+mytable.l.^2);
    case 'hexagonal'

    otherwise
        disp('Calculation is not available for this system.');
        ahkl=-1;
end

if ahkl~=-1
    [results,ahkl_T]=EvaluatePeakShift(sin2psi,ahkl,tth,hkl,GlobOpt);

    [~,ind]=sort(sin2psi);

    close all
    figure(1)
    hold on
    plot(sin2psi(ind),ahkl(ind),'bo','DisplayName','data');
    plot(sin2psi(ind),ahkl_T(ind),'r','DisplayName','fit');
    hold off
    xlabel('$$\sin^2\psi$$','interpreter','latex');
    ylabel('$$a_{hkl}$$ (\AA)','interpreter','latex');
    box on
    title([label ': Peak shift evaluation'])
    legend show
    set(gca,'Units','centimeters');
    pos=get(gca,'Position');
    set(gca,'Position',[pos(1:2) 12 8]);
    hL=findobj('type','legend');
    if results.sig_ef>0
        set(hL,'Location','southeast')
    else
        set(hL,'Location','northeast')
    end

    disp(results)
    DispResults(results,GlobOpt);
end

% switch PeakShiftModel
%     case 'sin2psi'
%         x0=[0 mean(ahkl)];
%         pT=lsqcurvefit(@(p,xdata) ahkl_rs(xdata,p,nu,E),x0,sin2psi,ahkl);
%         ahklT=ahkl_rs(sin2psi,pT,nu,E);
%         eqtext='$$a_{hkl}=a_0\left(\frac{\sigma}{E}\left[\left(1+\nu\right)\sin^2\psi-2\nu\right]+1\right)$$';
%         disp(['a0 = ' num2str(pT(1)) ' A']);
%         disp(['sig = ' num2str(pT(2)) ' GPa']);
%     case 'sin2psi_Theta'
%         x0=[0 mean(ahkl) 0];
%         pT=lsqcurvefit(@(p,xdata) ahkl_rs_samsh(xdata,p,nu,E,mytable.tth(ind)),x0,sin2psi,ahkl);
%         ahklT=ahkl_rs_samsh(sin2psi,pT,nu,E,mytable.tth(ind));
%         eqtext='$$a_{hkl}=a_0\left(\frac{\sigma}{E}\left[\left(1+\nu\right)\sin^2\psi-2\nu\right]+1\right)-a_0\cot\theta\textrm{d}\theta$$';
%         disp(['a0 = ' num2str(pT(1)) ' A']);
%         disp(['sig = ' num2str(pT(2)) ' GPa']);
%         disp(['dth = ' num2str(pT(3)*180/pi) ' deg']);
%     case 'VookWitt'
%         pT0=polyfit(sin2psi,ahkl,1);
%         x0=[pT0(2) 0 pT0(1) 0];
%         pT=lsqcurvefit(@(p,xdata) ahkl_VookWitt(xdata,p,[mytable.h(ind) mytable.k(ind) mytable.l(ind)]),x0,sin2psi,ahkl);
%         ahklT=ahkl_VookWitt(sin2psi,pT,[mytable.h(ind) mytable.k(ind) mytable.l(ind)]);
%         eqtext='$$a_{hkl}=P_1\sin^2\psi+P_2\Gamma\sin^2\psi+P_3\Gamma+P_4$$';
%     case 'Reuss'
%         pT0=polyfit(sin2psi,ahkl,1);
%         x0=[pT0(2) 0 0 pT0(1)];
%         pT=lsqcurvefit(@(p,xdata) ahkl_Reuss(xdata,p,[mytable.h(ind) mytable.k(ind) mytable.l(ind)],mytable.tth(ind)),x0,sin2psi,ahkl);
%         ahklT=ahkl_Reuss(sin2psi,pT,[mytable.h(ind) mytable.k(ind) mytable.l(ind)],mytable.tth(ind));
%         eqtext='$$a_{hkl}=P_1\sin^2\psi-3P_2\Gamma\sin^2\psi+2P_2\Gamma+P_3+P_4\cot\theta$$';
% end
% 
% [~,ind2]=sort(sin2psi);
% figure(1)
% hold on
% plot(sin2psi(ind2),ahkl(ind2),'bo','DisplayName','data');
% plot(sin2psi(ind2),ahklT(ind2),'r','DisplayName',eqtext);
% hold off
% rng=max(sin2psi)-min(sin2psi);
% for n=1:numel(ahkl)
%     text(sin2psi(ind2(n))+0.02*rng,ahkl(ind2(n)),[num2str(mytable.h(ind(ind2(n)))) ' ' num2str(mytable.k(ind(ind2(n)))) ' ' num2str(mytable.l(ind(ind2(n)))) ', grain ' mytable.grain{ind(ind2(n))}])
% end
% xlabel('$$\sin^2\psi$$','interpreter','latex');
% ylabel('$$a_{hkl}$$ (\AA)','interpreter','latex');
% box on
% title([label ': Peak shift evaluation'])
% legend show
% set(gca,'Units','centimeters');
% pos=get(gca,'Position');
% set(gca,'Position',[pos(1:2) 12 8]);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% PEAK BROADENING EVALUATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fwhm_meas=0.5*mytable.FWHM_left(ind)+0.5*mytable.FWHM_left(ind);
% fwhm_IF=sqrt((u0+u1*mytable.chi(ind)).*(tand(mytable.tth(ind)/2)).^2+(v0+v1*mytable.chi(ind)).*tand(mytable.tth(ind)/2)+(w0+w1*mytable.chi(ind)));
% % Gauss
% % betaG_meas=sqrt(pi/log(2))*fwhm_meas;
% % betaG_IF=sqrt(pi/log(2))*fwhm_IF;
% fwhmG=sqrt(fwhm_meas.^2-fwhm_IF.^2);
% betaG=1/lam*cosd(mytable.tth(ind)/2).*fwhmG*pi/360; % I need to recompute it to the rms of q-vector, fwhm is measured in 2th, therefore also devide it by 2, and measured in deg -> *pi/180, there is no 4*pi to be in correspondence with Scherrer equation
% % Cauchy
% fwhmC=fwhm_meas-fwhm_IF;
% betaC=1/lam*cosd(mytable.tth(ind)/2).*fwhmC*pi/360; % I need to recompute it to the rms of q-vector, fwhm is measured in 2th, therefore also devide it by 2, there is no 4*pi to be in correspondence with Scherrer equation
% 
% % beta=(1-mytable.eta(ind)).*betaG+mytable.eta(ind).*betaC;
% beta=betaG;
% 
% sinth=sind(mytable.tth(ind)/2);
% switch PeakBroadeningModel
%     case 'betahkl_classic'
%         x0=[10 0];
%         pT=lsqcurvefit(@(p,xdata) betahkl_classic(xdata,p,lam),x0,sinth,beta);
%         betahklT=betahkl_classic(sinth,pT,lam);
%         eqtext='$$\beta=\frac{0.94}{D}+\frac{4e}{\lambda}\sin\theta$$';
%         fprintf('D = %.3f nm\n',pT(1)/10);
%         fprintf('e = %.3e\n',pT(2));
%     case 'betahkl_disl_cubic'
%         x0=[10 0];
%         pT=lsqcurvefit(@(p,xdata) betahkl_disl_cubic(xdata,p,lam,C100,q,b,[mytable.h(ind) mytable.k(ind) mytable.l(ind)]),x0,sinth,beta);
%         betahklT=betahkl_disl_cubic(sinth,pT,lam,C100,q,b,[mytable.h(ind) mytable.k(ind) mytable.l(ind)]);
%         eqtext='$$\beta=\frac{0.94}{D}+\frac{4}{\lambda}\sqrt{2\pi C_{100}\left(1-q\Gamma\right)}bM\sqrt{\rho_{disl}}\sin\theta$$';
%         fprintf('D = %.3f nm\n',pT(1)/10);
%         fprintf('rho*M^2 = %.3e m^-2\n',(pT(2)^2)*1e20);
% end
% 
% [~,ind2]=sort(sinth);
% figure(2)
% hold on
% plot(sinth(ind2),beta(ind2),'bo','DisplayName','data');
% plot(sinth(ind2),betahklT(ind2),'r','DisplayName',eqtext);
% hold off
% rng=max(sinth)-min(sinth);
% for n=1:numel(sinth)
%     text(sinth(ind2(n))+0.02*rng,beta(ind2(n)),[num2str(mytable.h(ind(ind2(n)))) ' ' num2str(mytable.k(ind(ind2(n)))) ' ' num2str(mytable.l(ind(ind2(n)))) ', grain ' mytable.grain{ind(ind2(n))}])
% end
% xlabel('$$\sin\theta$$','interpreter','latex');
% ylabel('$$\beta$$ $$\left(\textrm{\AA}^{-1}\right)$$','interpreter','latex');
% box on
% title([label ': Williamson-Hall plot'])
% legend show
% set(gca,'Units','centimeters');
% pos=get(gca,'Position');
% set(gca,'Position',[pos(1:2) 12 8]);
% 
% hL=findobj('type','legend');
% set(hL,'Location','northwest','interpreter','latex');
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% AUXILIARY FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=ahkl_shift(x,p,hkl,tth,GlobOpt)
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
switch GlobOpt.PeakShiftModel
    case 'NoStress'
        % p(1)=a0
        y_stress=p(1);
    case 'Isotropic'
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
switch GlobOpt.sys
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

% shift influence (divergent beam and Bragg-Brentano geometry)
% p(6)=a0*dth
y_shift=p(6)./tand(tth/2);

y=y_stress+y_SF+y_shift;
end


function y=ahkl_ResultsEvaluation(pT,x,ahkl,hkl,GlobOpt)
h=hkl(:,1);
k=hkl(:,2);
l=hkl(:,3);
switch GlobOpt.PeakShiftModel
    case 'NoStress'
        % p(1)=a0
        y.a0=pT(1);
    case 'Isotropic'
        % a0*(sig*(1/2*s2*sin2psi+2*s1)+1)
        y.a_norm=pT(2);
        y.a_par=pT(1)+pT(2);
        y.sig_ef=pT(1); % a0*sig*(1+nu)/E
        if GlobOpt.Material.av==1
            y.sin2psi0=2*GlobOpt.Material.nu/(1+GlobOpt.Material.nu);
            y.a0=(y.a_par-y.a_norm)*y.sin2psi0+y.a_norm;
            y.sig=y.sig_ef*GlobOpt.Material.E/(y.a0.*(1+GlobOpt.Material.nu));
        end
    case 'VookWitt'
        % first the anisotropy needs to be substracted
        Gamma=(h.^2.*k.^2+k.^2.*l.^2+h.^2.*l.^2)./((h.^2+k.^2+l.^2).^2);
        ahkl_ani=pT(2).*Gamma.*x+pT(3).*Gamma;
        ahkl0=ahkl-ahkl_ani;
        pTlin=polyfit(x,ahkl0,1);
        y.a_norm=pTlin(2);
        y.a_par=pTlin(1)+pTlin(2);
        y.sig_ef=pTlin(1); % a0*sig*(1+nu_h00)/E_h00, but the stress needs to be evaluated from p(1) from the original fit, for cubic crystals E_h00=E and nu_h00=nu
        if GlobOpt.Material.av==1
            y.sin2psi0=2*GlobOpt.Material.nu/(1+GlobOpt.Material.nu);
            y.a0=(y.a_par-y.a_norm)*y.sin2psi0+y.a_norm;
            y.sig=pT(1)*GlobOpt.Material.E/(y.a0.*(1+GlobOpt.Material.nu));
        end
    case 'Reuss'
        % first the anisotropy needs to be substracted
        Gamma=(h.^2.*k.^2+k.^2.*l.^2+h.^2.*l.^2)./((h.^2+k.^2+l.^2).^2);
        ahkl_ani=-3*pT(2).*Gamma.*x+2*pT(2).*Gamma;
        ahkl0=ahkl-ahkl_ani;
        pTlin=polyfit(x,ahkl0,1);
        y.a_norm=pTlin(2);
        y.a_par=pTlin(1)+pTlin(2);
        y.sig_ef=pTlin(1); % a0*sig*(1+nu_h00)/E_h00, but the stress needs to be evaluated from p(1) from the original fit, for cubic crystals E_h00=E and nu_h00=nu
        if GlobOpt.Material.av==1
            y.sin2psi0=2*GlobOpt.Material.nu/(1+GlobOpt.Material.nu);
            y.a0=(y.a_par-y.a_norm)*y.sin2psi0+y.a_norm;
            y.sig=pT(1)*GlobOpt.Material.E/(y.a0.*(1+GlobOpt.Material.nu));
        end
end

if GlobOpt.FitSF==1
    % p(5)=a0*alpha
    y.alphaSF_ef=pT(5);
    if GlobOpt.Material.av==1
       y.alphaSF=pT(5)./y.a0;
    end
end

if GlobOpt.FitShift==1
    % p(5)=a0*dth
    y.dth_ef=pT(6);
    if GlobOpt.Material.av==1
       y.dth=pT(6)./y.a0;
    end
end
end

function [results,ahkl_T]=EvaluatePeakShift(sin2psi,ahkl,tth,hkl,GlobOpt)
p0=zeros(1,6);
fix=ones(1,6);
lb=-Inf*ones(1,6);
ub=+Inf*ones(1,6);
switch GlobOpt.PeakShiftModel
    case 'NoStress'
        p0(1)=mean(ahkl);
        fix(1:4)=[0 1 1 1];
        lb(1)=0.9*p0(1);
        ub(1)=1.1*p0(1);
    case 'Isotropic'
        pT0=polyfit(sin2psi,ahkl,1);
        p0(1:2)=pT0;
        fix(1:4)=[0 0 1 1];
        lb(1:2)=[0.8*p0(1) p0(2)/50];
        ub(1:2)=[1.2*p0(1) p0(2)*50];
    case 'VookWitt'
        pT0=polyfit(sin2psi,ahkl,1);
        p0(1:2)=pT0;
        fix(1:4)=[0 0 0 0];
    case 'Reuss'
        pT0=polyfit(sin2psi,ahkl,1);
        p0(1:2)=pT0;
        fix(1:4)=[0 0 0 1];
end

if GlobOpt.FitSF==1
    p0(5)=1e-2; % it is zero by default
    fix(5)=0;
    lb(5)=0;
    ub(5)=1; % is the alpha density or probability?
end

if GlobOpt.FitShift==1
    p0(6)=0;
    fix(6)=0;
    lb(6)=-0.1; % in radians
    ub(6)=0.1; % in radians
end

pT1=lsqcurvefit(@(p,xdata) ahkl_shift(xdata,interlace(p0,p,fix),hkl,tth,GlobOpt),p0(~fix),sin2psi,ahkl,lb(~fix),ub(~fix));
pT=interlace(p0,pT1,fix); % for the correct evaluation it is necessary to have all parameters, included those fixed to zero
results=ahkl_ResultsEvaluation(pT,sin2psi,ahkl,hkl,GlobOpt);
ahkl_T=ahkl_shift(sin2psi,pT,hkl,tth,GlobOpt);
end

function DispResults(results,GlobOpt)
if GlobOpt.Material.av==1
    fprintf(['a' char(8320) ' = %.4f ' char(8491) '\n'],results.a0);
    fprintf([char(963) ' = %.4f GPa\n'],results.sig);
else
    if strcmp(GlobOpt.PeakShiftModel,'NoStress')
        fprintf(['a' char(8320) ' = %.4f ' char(8491) '\n'],results.a0);
    else
        fprintf(['a' char(8741) ' = %.4f ' char(8491) '\n'],results.a_par);
        fprintf(['a' char(8869) ' = %.4f ' char(8491) '\n'],results.a_norm);
        fprintf(['a' char(8320) '*' char(963) '*(1+' char(957) ')/E = %.6f ' char(8491) '\n'],results.sig_ef);
    end    
end
end