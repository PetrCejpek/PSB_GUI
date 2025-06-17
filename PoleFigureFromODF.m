function PF=PoleFigureFromODF(hkl,ODF,sys)
% hkl - pole in reciprocal space
% ODF.TextureType - type of texture ('fibre',...')
% ODF.parameters - parameters of ODF teoretical formula corresponding to
%                  texture type specified above
% Examples:
% ODF.TextureType = 'fibre'
%   -> ODF.parameters.or1 - 1st orientation in angular (phi1,psi,phi2)
%                           space
%      ODF.parameters.or2 - 2nd orientation
%      ODF.parameters.width - width of texture function in degrees
% sys - crystal symmetry


% cubic
% I need to sum over the equivalent hkl indices
[V_equi,~]=GetEquivalentPlanes(hkl(:)',sys);

% prelocation
psi=0:3:90;
phi=0:5:360;
[PSI,PHI]=meshgrid(psi,phi);
PF=zeros(size(PSI));
for nv=1:size(V_equi)
    v=V_equi(nv,:)'/norm(hkl);
    % rotation corresponding to the center of the pole figure
    % - the diffraction vector must be along [0;0;1] in laboratory coordinates
    phi0=atan2d(v(2),v(1));
    psi0=acosd(v(3)./norm(v));
    R0=RotMat(phi0,psi0,0,'zxz');

    for m=1:size(PSI,1)
        for n=1:size(PSI,2)
            R=RotMat(0,PSI(m,n),PHI(m,n),'zxz');
            Q=R*R0*v; % in principle R0*v=[0;0;1]

            if sum(v==R0*v,3)
                phi1_m=ones(1,100)*PHI(m,n);
                psi_m=ones(1,100)*PSI(m,n);
                phi2_m=linspace(0,360,100);
            else
                lam=linspace(0,360,100);
                phi1_m=zeros(size(lam));
                psi_m=zeros(size(lam));
                phi2_m=zeros(size(lam));
                odfval=zeros(size(lam));
                for k=1:numel(lam)
                    RQ=RotMatArbitrary(Q,lam(k));
                    R2=RQ*R*R0;
 
                    % 'new' Euler angles (zxz) rotation
                    psi_m(k)=acosd(R2(3,3));
                    if R2(3,3)==1 | R2(3,3)==-1
                        phi1_m(k)=PHI(m,n);
                        phi2_m(k)=mod(atan2d(R2(2,1),R2(1,1)),360)-PHI(m,n);
                    elseif R2(3,3)==-1
                        phi1_m(k)=PHI(m,n);
                        phi2_m(k)=mod(atan2d(R2(2,1),R2(1,1)),360)+PHI(m,n);
                    else
                        phi2_m(k)=mod(atan2d(R2(1,3)/sind(psi_m(k)),-R2(2,3)/sind(psi_m(k))),360);
                        phi1_m(k)=mod(atan2d(R2(3,1)/sind(psi_m(k)),R2(3,2)/sind(psi_m(k))),360);
                    end
                end
            end
            odfval=ODFevaluate(phi1_m,psi_m,phi2_m,ODF,sys);
            PF(m,n)=PF(m,n)+sum(odfval);

        % plot of the measurement in 3D Euler space (phi1 around z, psi around x'',
        % phi2 around z'')
%         plot3(phi1_m,psi_m,phi2_m)
%         axis equal
%         xlim([0 360])
%         ylim([0 180])
%         zlim([0 360])
%         box on
%         xlabel('\Phi_1 (deg)')
%         ylabel('\Psi (deg)')
%         zlabel('\Phi_2 (deg)')
%         title(['hkl = ' num2str(hkl(1)) ' ' num2str(hkl(2)) ' ' num2str(hkl(3)) ', \Psi = ' num2str(PSI(m,n)) ' deg, \Phi = ' num2str(PHI(m,n)) ' deg'])
%         input('')
        end
    end
end

pcolor(PSI.*cosd(PHI),PSI.*sind(PHI),PF)
shading flat