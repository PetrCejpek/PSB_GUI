function [E,nu,G]=ElasticConstantsInDirection(Sij,m1,m2,m3)
% I have 'testing' uniaxial stress with the magnitude 1 along direction m
% m could be a list (n x 3)
nm=sqrt(m1.*m1+m2.*m2+m3.*m3);
m1=m1./nm; m2=m2./nm; m3=m3./nm;

E=zeros(size(m1));
nu=zeros(size(m1));
G=zeros(size(m1));
for n=1:numel(m1)
    sig=[m1(n).*m1(n);m2(n).*m2(n);m3(n).*m3(n);0.5*m2(n).*m3(n);0.5*m1(n).*m3(n);0.5*m1(n).*m2(n)];
    eps=Sij*sig;

    % Now I have to find the strain in the correct direction
    % For Young modulus, the 'testing' direction of strain should be the same
    % as for stress
    % 1/E=eps_m/norm(sig), but norm(sig)=1
    eps_m=m1(n).*m1(n).*eps(1)+...
          m2(n).*m2(n).*eps(2)+...
          m3(n).*m3(n).*eps(3)+...
          2*m2(n).*m3(n).*eps(4)+...
          2*m1(n).*m3(n).*eps(5)+...
          2*m1(n).*m2(n).*eps(6);
    E(n)=(1./eps_m);

    % For nu, we have 2 'testing' directions of strain perpendicular to the
    % direction of stress
    if m1(n)==0 && m2(n)==0
        v1=[1;0;0];
        v2=[0;1;0];
    elseif m1(n)==0 && m3(n)==0
        v1=[1;0;0];
        v2=[0;0;1];
    elseif m2(n)==0 && m3(n)==0
        v1=[0;1;0];
        v2=[0;0;1];
    elseif m3(n)==0
        v1=[1;-m1(n)/m2(n),1];
        v1=v1./norm(v1);
        RQ=[m1(n)^2+(1-m1(n)^2)*cosd(90), m1(n)*m2(n)*(1-cosd(90))-m3(n)*sind(90), m1(n)*m3(n)*(1-cosd(90))+m2(n)*sind(90);
            m1(n)*m2(n)*(1-cosd(90))+m3(n)*sind(90), m2(n)^2+(1-m2(n)^2)*cosd(90), m2(n)*m3(n)*(1-cosd(90))-m1(n)*sind(90);
            m1(n)*m3(n)*(1-cosd(90))-m2(n)*sind(90), m2(n)*m3(n)*(1-cosd(90))+m1(n)*sind(90), m3(n)^2+(1-m3(n)^2)*cosd(90)];
        v2=RQ*v1;
    else
        v1=[1;1;-(m1(n)+m2(n))/m3(n)];
        v1=v1./norm(v1);
        % rotation of 90 deg around the measurement direction
        RQ=[m1(n)^2+(1-m1(n)^2)*cosd(90), m1(n)*m2(n)*(1-cosd(90))-m3(n)*sind(90), m1(n)*m3(n)*(1-cosd(90))+m2(n)*sind(90);
            m1(n)*m2(n)*(1-cosd(90))+m3(n)*sind(90), m2(n)^2+(1-m2(n)^2)*cosd(90), m2(n)*m3(n)*(1-cosd(90))-m1(n)*sind(90);
            m1(n)*m3(n)*(1-cosd(90))-m2(n)*sind(90), m2(n)*m3(n)*(1-cosd(90))+m1(n)*sind(90), m3(n)^2+(1-m3(n)^2)*cosd(90)];
        v2=RQ*v1;
    end
    eps_v1=v1(1).*v1(1).*eps(1)+...
           v1(2).*v1(2).*eps(2)+...
           v1(3).*v1(3).*eps(3)+...
           2*v1(2).*v1(3).*eps(4)+...
           2*v1(1).*v1(3).*eps(5)+...
           2*v1(1).*v1(2).*eps(6);
    eps_v2=v2(1).*v2(1).*eps(1)+...
           v2(2).*v2(2).*eps(2)+...
           v2(3).*v2(3).*eps(3)+...
           2*v2(2).*v2(3).*eps(4)+...
           2*v2(1).*v2(3).*eps(5)+...
           2*v2(1).*v2(2).*eps(6);
    nu(n)=-(eps_v1+eps_v2)/(2*eps_m);

    % G=gamma_xy/tau_xy, which in this case should correspond to eps(4:6)
    % and sig(4:6)
    % New 'testing' stress
%     sig=[0;0;0;m2(n)*m3(n);m1(n)*m3(n);m1(n)*m2(n)];
%     eps=Sij*sig;

    % Now I have to find the strain in the correct direction
    % For Young modulus, the 'testing' direction of strain should be the same
    % as for stress
    % 1/E=eps_m/norm(sig), but norm(sig)=1
%     eps_m=2*m2(n).*m3(n).*eps(4)+...
%           2*m1(n).*m3(n).*eps(5)+...
%           2*m1(n).*m2(n).*eps(6);
    G(n)=sum(sig(4:6))/sum(eps(4:6));
end