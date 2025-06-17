function C=GetFullElasticTensor(Cr)
% C is the 3x3x3x3 tensor/matrix
% Cr is the 6x6 tensor/matrix

%C_ijkl=C_jikl
%C_ijkl=C_ijlk

% Cr=[C(1,1,1,1) C(1,1,2,2) C(1,1,3,3) C(1,1,2,3) C(1,1,1,3) C(1,1,1,2);
%     C(2,2,1,1) C(2,2,2,2) C(2,2,3,3) C(2,2,2,3) C(2,2,1,3) C(2,2,1,2);
%     C(3,3,1,1) C(3,3,2,2) C(3,3,3,3) C(3,3,2,3) C(3,3,1,3) C(3,3,1,2);
%     C(2,3,1,1) C(2,3,2,2) C(2,3,3,3) C(1,2,2,3) C(2,3,1,3) C(2,3,1,2);
%     C(1,3,1,1) C(1,3,2,2) C(1,3,3,3) C(1,3,2,3) C(1,3,1,3) C(1,3,1,2);
%     C(1,2,1,1) C(1,2,2,2) C(1,2,3,3) C(2,3,2,3) C(1,2,1,3) C(1,2,1,2);];

C=zeros(3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3

                if i==j
                    m=i;
                else
                    m=9-(i+j); % 1 2, 2 1 -> 6; 1 3, 3 1 -> 5; 2 3, 3 2 -> 4
                end
                
                if k==l
                    n=k;
                else
                    n=9-(k+l); % 1 2, 2 1 -> 6; 1 3, 3 1 -> 5; 2 3, 3 2 -> 4
                end

                C(i,j,k,l)=Cr(m,n);
            end
        end
    end
end