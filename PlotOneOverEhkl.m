function PlotOneOverEhkl(Sij)
% Sij is 6x6 compliance matrix

[PSI,PHI]=meshgrid(linspace(0,180,60),linspace(0,360,72)); % grid with 3 deg in psi and 5 deg in phi
VAL=OneOverEhkl(Sij,PSI,PHI);

X=1./VAL.*sind(PSI).*cosd(PHI);
Y=1./VAL.*sind(PSI).*sind(PHI);
Z=1./VAL.*cosd(PSI);

surf(X,Y,Z,'CData',1./VAL);
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
title('$$E_{hkl}$$','interpreter','latex');
colorbar