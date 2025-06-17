function [phi1,psi,phi2]=GetEulerAnglesZXZ(R,rotopt)
% R - rotation/transformation matrix
% rotopt - 'extrinsic' for extrinsic rotation
%        - 'intrinsic' for intrinsic rotation

psi=acosd(R(3,3));
if R(3,3)==1 | R(3,3)==-1
    phi1=0;
    phi2=mod(atan2d(R(2,1),R(1,1)),360);
elseif R(3,3)==-1
    phi1=0;
    phi2=mod(atan2d(R(2,1),R(1,1)),360);
else
    switch rotopt
        case 'extrinsic'
            % R=R3*R2*R1
            phi2=mod(atan2d(R(1,3)/sind(psi),-R(2,3)/sind(psi)),360);
            phi1=mod(atan2d(R(3,1)/sind(psi),R(3,2)/sind(psi)),360);
        case 'intrinsic'
            % R=R1*R2*R3
            phi1=mod(atan2d(R(1,3)/sind(psi),-R(2,3)/sind(psi)),360);
            phi2=mod(atan2d(R(3,1)/sind(psi),R(3,2)/sind(psi)),360);
    end
end