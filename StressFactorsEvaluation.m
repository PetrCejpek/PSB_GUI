function F=StressFactorsEvaluation(s_ij,c_ij,ODFfilename,method,hkl,psim,phim,geom,sys,ac,plotopt)
% s_ij compliance matrix (6x6)
% c_ij stiffness matrix (6x6)
% ODFfilename - name of the file containing ODF data in scattered format (4 columns)
% method - name of the method to compute stress factors ('Reuss', 'Modified
%          Voigt', 'Kroener', 'Inverse Kroener')
% hkl - diffraction indices of the measured data
% psim - psi angles of the measured data
% phim - phi angles of the measured data
% geom - measurement geometry ('GAXRD', 'Pole figure')
% sys - crystallographic system ('cubic', 'hexagonal')
% ac - lattice parameters (for cubic system can be ignored, for hexagonal a=ac(1), c=ac(2)
% plotopt - if 1, "integration trails" through the Euler space are plotted

N=size(hkl,1); % number of measured Bragg peaks

%% Loading of ODF data if available
if ODFfilename~=1
    % if ODFfilename is 1, the case of ideal polycrystal is taken into
    % acount (ODF=1 for all of angles)
    ODFdata=importdata(ODFfilename); % data should be in 4 columns -> Phi1, Psi, Phi2, ODF
    ODFinterpolant=scatteredInterpolant(ODFdata.data(:,1),ODFdata.data(:,2),ODFdata.data(:,3),ODFdata.data(:,4));
end

%% Computation of angles related to hkl
switch sys
    case 'cubic'
        PsiB=acosd(hkl(:,3)./sqrt(hkl(:,1).^2+hkl(:,2).^2+hkl(:,3).^2)); % angle between hkl and direction [001]
%         betaB=acosd(hkl(:,1)./sqrt(hkl(:,1).^2+hkl(:,2).^2));
        betaB=atan2d(hkl(:,2),hkl(:,1)); % angle between projection of hkl to (001) and direction [100]
    case 'hexagonal'
        hx=hkl(:,1)/ac(1);
        hy=-hkl(:,1)./ac(1).*tand(120)+hkl(:,2)./(ac(1).*sind(120));
        hz=hkl(:,3)./ac(2);
        PsiB=acosd(hz./sqrt(hx.^2+hy.^2+hz.^2));
%         betaB=acosd(hx./sqrt(hx.^2+hy.^2));
        betaB=atan2d(hy,hx);
    otherwise
        disp('Uknown system!')
end

%% Computation of the trajectory for integration in Euler space
lam=linspace(0,360,361);
Phi1=zeros(N,numel(lam));
Psi=zeros(N,numel(lam));
Phi2=zeros(N,numel(lam));
for m=1:N
    % one need to compute the orientation of crystallite, when the Qhkl
    % vector is pointing to the detector
    % G1 - rotation matrix, which rotate the crystal to the position, when
    %      Qhkl is pointing along [001] global/laboratory axis, then it
    %      does additional rotation along this direction an angle lam
    % G2 - rotates the crystal to the position, when Qhkl is pointing to
    %      the detector
    %______________________________________________________________________
    % the second rotation, rotations are extrinsic
    switch geom
        case 'GAXRD'
            G2=RotMat(0,psim(m),phim(m),'zyz'); % for GAXRD setup
        case 'Pole figure'
            G2=RotMat(0,psim(m),phim(m),'zxz'); % for Pole figure-like setup
    end
    for n=1:numel(lam)
        % the first rotation, rotations are extrinsic
        % FOR CONTROL: if we have a single-crystal and we are measuring
        % directly at hkl position in pole-figure like geometry, phim
        % should compensate with betaB and psim with PsiB 
        G1=RotMat(-betaB(m),-PsiB(m),lam(n),'zyz');
        R=G2*G1;
        [Phi1(m,n),Psi(m,n),Phi2(m,n)]=GetEulerAnglesZXZ(R,'intrinsic'); % IMPORTANT: these angles should be intrinsic!
    end
end

%% Plot trails
if plotopt==1
    figure('name','ODF trails')
    clrs='bgrcmyk';
    hold on
    for m=1:size(Phi1,1)
        plot3(Phi1(m,:),Psi(m,:),Phi2(m,:),'Color',clrs(mod(m-1,7)+1),'DisplayName',['Measurement #' num2str(m)])
    end
    x0=min(ODFdata.data(:,1));
    x1=max(ODFdata.data(:,1));
    y0=min(ODFdata.data(:,2));
    y1=max(ODFdata.data(:,2));
    z0=min(ODFdata.data(:,3));
    z1=max(ODFdata.data(:,3));
    plot3([x0 x1 x1 x0 x0 x0 x1 x1 x0 x0 x1 x1 x1 x1 x0 x0 x0],...
          [y0 y0 y1 y1 y0 y0 y0 y1 y1 y0 y0 y0 y1 y1 y1 y1 y0],...
          [z0 z0 z0 z0 z0 z1 z1 z1 z1 z1 z0 z1 z0 z1 z0 z1 z0],...
          'k--','DisplayName','Limits of ODF values')
    hold off
    axis equal
    xlim([0 360])
    ylim([0 180])
    zlim([0 360])
    xlabel('\Phi_1 (deg)')
    ylabel('\Psi (deg)')
    zlabel('\Phi_2 (deg)')
    legend show
end

%% ODF evaluation
odfval=zeros(N,numel(lam));
for m=1:N
    if ODFfilename==1
        % ideal polycrystal
        odfval(m,:)=ones(size(lam));
    else
        switch sys
            case 'cubic'
                odfval(m,:)=ODFinterpolant(Phi1(m,:),90-abs(Psi(m,:)-90),mod(Phi2(m,:),90));
                %         odfval(m,:)=ODFinterpolant(Phi1(m,:),Psi(m,:),Phi2(m,:));
            case 'hexagonal'
                % IT NEEDS TO BE CHECKED
                odfval(m,:)=ODFinterpolant(Phi1(m,:),90-abs(Psi(m,:)-90),mod(Phi2(m,:),120));
        end
    end
end

%% Own model computation
F=zeros(3,3,size(hkl,1));
Cfull=GetFullElasticTensor(c_ij);
Sfull=GetFullComplianceTensor(s_ij);
switch method
    case 'Kroener'

    case 'InverseKroener'

    case 'ModifiedVoigt'
        A=zeros(N,1);
        for m=1:N
            A(m)=TrapezoidIntegral(lam,odfval(m,:)); % normalisation
            Crot=zeros(3,3,3,3,numel(lam)); % all rotated stiffness tensors
%             Srot=zeros(3,3,3,3,numel(lam)); % all rotated compliance tensors
            measd=[sind(psim(m))*cosd(phim(m)); sind(psim(m))*sind(phim(m)); cosd(psim(m))];
            % computation of every rotated full elastic tensor beforehand
            for n=1:numel(lam)
                R=RotMat(Phi2(m,n),Psi(m,n),Phi1(m,n),'zxz'); % because we have intrinsic rotation, I shoul switch Phi1 and Phi2
                Crot(:,:,:,:,n)=RotateFullElasticTensor(Cfull,R); % one of this may not be neccessary depending on method
%                 Srot(:,:,:,:,n)=RotateFullElasticTensor(Sfull,R); % one of this may not be neccessary depending on method
            end

            Cav=zeros(3,3,3,3,numel(lam));
            for i=1:3
                for j=1:3
                    % let's compute the function inside the integral
                    y=zeros(size(lam));
                    for k=1:3
                        for l=1:3
                            hlp=Crot(i,j,k,l,:);
                            y=hlp(:)';

                            % integrate and compute stress factors
                            Cav(i,j,k,l,m)=TrapezoidIntegral(lam,y.*odfval(m,:))/A(m);
                        end
                    end
                end
            end
            % get inverse of the average tensor
            c=GetReducedElasticTensor(Cav);
            ci=inv(c);
            Ci=GetFullStiffnessTensor(ci);
            for i=1:3
                for j=1:3
                    for k=1:3
                        for l=1:3
                            F(i,j)=F(i,j)+Ci(i,j,k,l)*measd(k)*measd(l);
                        end
                    end
                end
            end
        end
    case 'Reuss'
        A=zeros(N,1);
        for m=1:N
            A(m)=TrapezoidIntegral(lam,odfval(m,:)); % normalisation
%             Crot=zeros(3,3,3,3,numel(lam)); % all rotated stiffness tensors
            Srot=zeros(3,3,3,3,numel(lam)); % all rotated compliance tensors
            measd=[sind(psim(m))*cosd(phim(m)); sind(psim(m))*sind(phim(m)); cosd(psim(m))];
            % computation of every rotated full elastic tensor beforehand
            for n=1:numel(lam)
                R=RotMat(Phi2(m,n),Psi(m,n),Phi1(m,n),'zxz'); % because we have intrinsic rotation, I shoul switch Phi1 and Phi2
%                 Crot(:,:,:,:,n)=RotateFullElasticTensor(Cfull,R); % one of this may not be neccessary depending on method
                Srot(:,:,:,:,n)=RotateFullElasticTensor(Sfull,R); % one of this may not be neccessary depending on method
            end

            for i=1:3
                for j=1:3
                    % let's compute the function inside the integral
                    y=zeros(size(lam));
                    for k=1:3
                        for l=1:3
                            hlp=Srot(i,j,k,l,:)*measd(k)*measd(l);
                            y=y+hlp(:)';
                        end
                    end

                    % integrate and compute stress factors
                    F(i,j,m)=TrapezoidIntegral(lam,y.*odfval(m,:))/A(m);
                end
            end
        end
end