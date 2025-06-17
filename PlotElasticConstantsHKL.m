function PlotElasticConstantsHKL(Sij,R)
% Sij is 6x6 compliance matrix
% R is rotation matrix 3x3

[PSI,PHI]=meshgrid(linspace(0,180,60),linspace(0,360,72)); % grid with 3 deg in psi and 5 deg in phi
H=sind(PSI).*cosd(PHI);
K=sind(PSI).*sind(PHI);
L=cosd(PSI);
if nargin>1
    Sij1=GetReducedComplianceTensor(RotateFullElasticTensor(GetFullComplianceTensor(Sij),R));
    Sij2=inv(GetReducedElasticTensor(RotateFullElasticTensor(GetFullElasticTensor(inv(Sij)),R)));
    sum(abs(Sij1-Sij2),'all')
    Sij=Sij2;
%     Sij=inv(RotateStiffnessTensor(R,inv(Sij)));
%     for n=1:numel(H)
%         vec=R*[H(n);K(n);L(n)];
%         H(n)=vec(1);
%         K(n)=vec(2);
%         L(n)=vec(3);
%     end
end


% [E,NU,G]=ElasticConstantsHKL(Sij,H,K,L);
[E,NU,G]=ElasticConstantsInDirection(Sij,H,K,L);

figure('Name','E_hkl')
X=E.*sind(PSI).*cosd(PHI);
Y=E.*sind(PSI).*sind(PHI);
Z=E.*cosd(PSI);
max(E,[],'all')
min(E,[],'all')
surf(X,Y,Z,'CData',E);
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
title('$$E_{hkl}$$','interpreter','latex');
colorbar

figure('Name','nu_hkl')
X=NU.*sind(PSI).*cosd(PHI);
Y=NU.*sind(PSI).*sind(PHI);
Z=NU.*cosd(PSI);

surf(X,Y,Z,'CData',NU);
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
title('$$\nu_{hkl}$$','interpreter','latex');
colorbar

figure('Name','G_hkl')
X=G.*sind(PSI).*cosd(PHI);
Y=G.*sind(PSI).*sind(PHI);
Z=G.*cosd(PSI);

surf(X,Y,Z,'CData',G);
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
title('$$G_{hkl}$$','interpreter','latex');
colorbar