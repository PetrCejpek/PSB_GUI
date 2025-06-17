function [planes,m_hkl]=GetEquivalentPlanes(hkl,sys)

switch sys
    case 'cubic'
        I=perms(1:3);

        % hkl
        HKL=hkl(I);

        % -h k l
        hkl1=[-hkl(1) hkl(2) hkl(3)];
        HKL=[HKL; hkl1(I)];

        % h -k l
        hkl1=[hkl(1) -hkl(2) hkl(3)];
        HKL=[HKL; hkl1(I)];

        % h k -l
        hkl1=[hkl(1) hkl(2) -hkl(3)];
        HKL=[HKL; hkl1(I)];

        % -h -k l
        hkl1=[-hkl(1) -hkl(2) hkl(3)];
        HKL=[HKL; hkl1(I)];

        % -h k -l
        hkl1=[-hkl(1) hkl(2) -hkl(3)];
        HKL=[HKL; hkl1(I)];

        % h -k -l
        hkl1=[hkl(1) -hkl(2) -hkl(3)];
        HKL=[HKL; hkl1(I)];

        % -h -k -l
        hkl1=[-hkl(1) -hkl(2) -hkl(3)];
        HKL=[HKL; hkl1(I)];

        planes=unique(HKL,'rows');
        m_hkl=size(planes,1);
    case 'hexagonal'

    otherwise
        disp('Uknown option.')
end