function [ODF,PHI,PSI,LAM]=FibreTextureODF(or1,or2,sigma)
tic
% coordinate limits
phi=0:5:360;
psi=0:5:90;
lam=0:5:90;

% coordinates
[PHI,PSI,LAM]=meshgrid(phi,psi,lam);

% distance from the line connecting or1 and or2
d=zeros(size(PHI));
for n=1:numel(PHI)
    d(n)=norm(cross(or1-[PHI(n) PSI(n) LAM(n)],or2-[PHI(n) PSI(n) LAM(n)]))/norm(or1-or2);
end

F=exp(-d.^2/(2*sigma)^2);

ODF=F/sum(F/(8*pi^2),'all');
toc

tic
ind=(d<=sigma);
plot3(PHI(ind),PSI(ind),LAM(ind),'o')
axis equal
xlim([0 360]);
ylim([0 90]);
zlim([0 90]);
xlabel('\Phi (deg)');
ylabel('\Psi (deg)');
zlabel('\Lambda (deg)');
box
toc