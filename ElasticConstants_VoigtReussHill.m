function [nu,E,G,B]=ElasticConstants_VoigtReussHill(C,system)
% Computation of mechanical properties nu, E, G, B from stiffness tensor
% with Voigt-Reuss-Hill model
switch system
    case {'cubic','cubic_fcc','cubic_bcc'}
        BV=(C(1,1)+2*C(1,2))/3;
        GV=(C(1,1)-C(1,2)+3*C(4,4))/5;
        BR=BV;
        GR=5*(C(1,1)-C(1,2))*C(4,4)/(3*C(1,1)-3*C(1,2)+4*C(4,4));
    case 'hexagonal'
        BV=2*(C(1,1)+C(1,2)+C(3,3)/2+2*C(1,3))/9;
        GV=(7*C(1,1)-5*C(1,2)+12*C(4,4)+2*C(3,3)-4*C(1,3))/30;

        M=C(1,1)+C(1,2)+2*C(3,3)-4*C(1,3);
        C2=(C(1,1)+C(1,2))*C(3,3)-2*C(1,3)^2;
        BR=C2/M;
        GR=5/2*(C2*C(4,4)*C(6,6)/(3*BV*C(4,4)*C(6,6)+C2*(C(4,4)+C(6,6))));
end

B=(BR+BV)/2;
G=(GR+GV)/2;
E=9*B*G/(3*B+G);
nu=(3*B-2*G)/(6*B+2*G);