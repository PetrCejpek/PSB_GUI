function R=RotMat(al,be,ga,opt,opt2)
% The function will compute the rotation matrix with the 3 al, be and ga in
% degrees.
% The axes around which the rotation happens are set with the opt parameter
% as the triplet of x, y and z.
% Example: RotMat(al,be,ga,'zxz') is the rotation matrix with classical
% Eulerian angles al, be, ga. The first rotation is around z with angle al,
% the second rotation around x' with angle be, and the third rotation is
% around z'' with angle ga.
% opt2 says, if the rotation is extrinsic o intrisic
if nargin<5
    % default rotation type
    opt2='extrinsic';
end

if (contains(opt,'x') | contains(opt,'y') | contains(opt,'z')) & length(opt)==3
    if opt(1)==opt(2) | opt(2)==opt(3)
        warning('Consecutive rotations should not be around the same axis.');
    end

    switch opt(1)
        case 'x'
            R1=[1 0 0; 0 cosd(al) -sind(al); 0 sind(al) cosd(al)];
        case 'y'
            R1=[cosd(al) 0 -sind(al); 0 1 0; sind(al) 0 cosd(al)];
        case 'z'
            R1=[cosd(al) -sind(al) 0; sind(al) cosd(al) 0; 0 0 1];
    end

    switch opt(2)
        case 'x'
            R2=[1 0 0; 0 cosd(be) -sind(be); 0 sind(be) cosd(be)];
        case 'y'
            R2=[cosd(be) 0 -sind(be); 0 1 0; sind(be) 0 cosd(be)];
        case 'z'
            R2=[cosd(be) -sind(be) 0; sind(be) cosd(be) 0; 0 0 1];
    end

    switch opt(3)
        case 'x'
            R3=[1 0 0; 0 cosd(ga) -sind(ga); 0 sind(ga) cosd(ga)];
        case 'y'
            R3=[cosd(ga) 0 -sind(ga); 0 1 0; sind(ga) 0 cosd(ga)];
        case 'z'
            R3=[cosd(ga) -sind(ga) 0; sind(ga) cosd(ga) 0; 0 0 1];
    end

    switch opt2
        case 'extrinsic'
            R=R3*R2*R1;
        case 'intrinsic'
            R=R1*R2*R3;
        otherwise
            warning('Uknown rotation type.');
            R=NaN;
    end
else
    warning('Uknown axes option.');
    R=NaN;
end