function T=UnfoldComplianceTensor_2_to_4(C)
ind=[1 1;
    2 2;
    3 3;
    2 3;
    1 3;
    1 2;];
% possible permutations
per=[1 2 3 4; % ijkl
    2 1 3 4; % jikl
    2 1 4 3; % jilk
    1 2 4 3; % ijlk
    3 4 1 2; % klij
    4 3 1 2; % lkij
    4 3 2 1; % lkji
    3 4 2 1]; % klji

T=zeros(3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                hlp=[i j k l];
                for x=1:size(per,1)
                    % choose permutation from per
                    i1=hlp(per(x,1));
                    j1=hlp(per(x,2));
                    k1=hlp(per(x,3));
                    l1=hlp(per(x,4));
                    % find values in mapping indices
                    p=find(i1==ind(:,1) & j1==ind(:,2));
                    q=find(k1==ind(:,1) & l1==ind(:,2));
                    if ~isempty(p) && ~isempty(q)
                        break;
                    end
                end
                T(i,j,k,l)=C(p,q);
            end
        end
    end
end