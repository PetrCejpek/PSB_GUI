% R1=RotMat(10,20,30,'zxz');
% 
% phi=0;
% psi=20;
% lam=30;
% 
% m=[sind(psi).*cosd(phi); sind(psi).*sind(phi); cosd(psi)];
% RQ=[m(1)^2+(1-m(1)^2)*cosd(lam), m(1)*m(2)*(1-cosd(lam))-m(3)*sind(lam), m(1)*m(3)*(1-cosd(lam))+m(2)*sind(lam);
%     m(1)*m(2)*(1-cosd(lam))+m(3)*sind(lam), m(2)^2+(1-m(2)^2)*cosd(lam), m(2)*m(3)*(1-cosd(lam))-m(1)*sind(lam);
%     m(1)*m(3)*(1-cosd(lam))-m(2)*sind(lam), m(2)*m(3)*(1-cosd(lam))+m(1)*sind(lam), m(3)^2+(1-m(3)^2)*cosd(lam)];
% 
% R2=RQ*R1;
% 
% PSI=acosd(R2(3,3));
% PHI2=atan2d(R2(1,3)/sind(PSI),-R2(2,3)/sind(PSI));
% PHI1=atan2d(R2(3,1)/sind(PSI),R2(3,2)/sind(PSI));
% 
% Rnew=RotMat(PHI1,PSI,PHI2,'zxz');
% Rnew-R2

phi=0; % measurement position, or position in PF
psi=0; % measurement position, or position in PF
lam=linspace(0,359,100); % angles for rotation around Q
m=[sind(psi).*cosd(phi); sind(psi).*sind(phi); cosd(psi)]; % measurement direction, or direction of Qhkl

phi1_m=zeros(size(lam));
psi_m=zeros(size(lam));
phi2_m=zeros(size(lam));
for n=1:numel(lam)
    % initial rotation of crystallite axes with respect to laboratory
    % system
    R1=RotMat(phi,psi,0,'zxz');

    % additional rotation around measurement direction
    RQ=[m(1)^2+(1-m(1)^2)*cosd(lam(n)), m(1)*m(2)*(1-cosd(lam(n)))-m(3)*sind(lam(n)), m(1)*m(3)*(1-cosd(lam(n)))+m(2)*sind(lam(n));
        m(1)*m(2)*(1-cosd(lam(n)))+m(3)*sind(lam(n)), m(2)^2+(1-m(2)^2)*cosd(lam(n)), m(2)*m(3)*(1-cosd(lam(n)))-m(1)*sind(lam(n));
        m(1)*m(3)*(1-cosd(lam(n)))-m(2)*sind(lam(n)), m(2)*m(3)*(1-cosd(lam(n)))+m(1)*sind(lam(n)), m(3)^2+(1-m(3)^2)*cosd(lam(n))];

    % resulting rotation matrix
    R2=RQ*R1;

    % 'new' Euler angles (zxz) rotation
    psi_m(n)=acosd(R2(3,3));
    if R2(3,3)==1 | R2(3,3)==-1
        phi1_m(n)=phi;
        phi2_m(n)=mod(atan2d(R2(2,1),R2(1,1)),360)-phi;
    elseif R2(3,3)==-1
        phi1_m(n)=phi;
        phi2_m(n)=mod(atan2d(R2(2,1),R2(1,1)),360)+phi;
    else
        phi2_m(n)=mod(atan2d(R2(1,3)/sind(psi_m(n)),-R2(2,3)/sind(psi_m(n))),360);
        phi1_m(n)=mod(atan2d(R2(3,1)/sind(psi_m(n)),R2(3,2)/sind(psi_m(n))),360);
    end
end

% plot of the measurement in 3D Euler space (phi1 around z, psi around x'',
% phi2 around z'')
plot3(phi1_m,psi_m,phi2_m)
axis equal
xlim([0 360])
ylim([0 180])
zlim([0 360])
box on
xlabel('\Phi_1 (deg)')
ylabel('\Psi (deg)')
zlabel('\Phi_2 (deg)')