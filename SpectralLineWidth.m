function w=SpectralLineWidth(tth,U,V,W,lam)
% Function w=SpectralLineWidth(tth,U,V,W) will return the width of spectral
% line (in Angstroems) for the Bragg angle tth (in deg) and an instrumental
% broadening given by the parameters U, V, W of Caglioti polynomial (also
% in deg).
% tth - Bragg angle (deg)
% U, V, W - Caglioti polynomial parameters (deg)
% lam - wavelength of the used X-ray radiation (Angstroem)

% instrumental broadening in deg
FWHM=sqrt(U.*tand(tth/2).^2+V.*(tand(tth/2))+W);

% recompute into Angstroem
% FWHM=lam/(4*pi).*cosd(tth/2)./(sind(tth/2).^2).*FWHM*pi/360;
FWHM=1/lam.*cosd(tth/2).*FWHM*pi/360;

w=lam.^2./FWHM;
% w=1./FWHM;