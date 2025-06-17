function S=GetFullComplianceTensor(Sr)
% S is the 3x3x3x3 tensor/matrix
% Sr is the 6x6 tensor/matrix

%C_ijkl=C_jikl
%C_ijkl=C_ijlk
%S_ijkl=S_jikl
%S_ijkl=S_ijlk
% S_ijkl=S_mn % m,n are 1,2,3
% S_ijkl=S_mn/2 % m is 1,2,3, n is 4,5,6 or vice versa
% S_ijkl=S_mn/4 % both m,n are 4,5,6

% Sr=[S(1,1,1,1) S(1,1,2,2) S(1,1,3,3) S(1,1,1,2) S(1,1,1,3) S(1,1,2,3);
%     S(2,2,1,1) S(2,2,2,2) S(2,2,3,3) S(2,2,1,2) S(2,2,1,3) S(2,2,2,3);
%     S(3,3,1,1) S(3,3,2,2) S(3,3,3,3) S(3,3,1,2) S(3,3,1,3) S(3,3,2,3);
%     S(1,2,1,1) S(1,2,2,2) S(1,2,3,3) S(1,2,1,2) S(1,2,1,3) S(1,2,2,3);
%     S(1,3,1,1) S(1,3,2,2) S(1,3,3,3) S(1,3,1,2) S(1,3,1,3) S(1,3,2,3);
%     S(2,3,1,1) S(2,3,2,2) S(2,3,3,3) S(2,3,1,2) S(2,3,1,3) S(2,3,2,3);];

S=zeros(3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3

                if i==j
                    m=i; % 1 1 -> 1; 2 2 -> 2; 3 3 -> 3
                else
                    m=9-(i+j); % 1 2, 2 1 -> 6; 1 3, 3 1 -> 5; 2 3, 3 2 -> 4
                end
                
                if k==l
                    n=k; % 1 1 -> 1; 2 2 -> 2; 3 3 -> 3
                else
                    n=9-(k+l); % 1 2, 2 1 -> 6; 1 3, 3 1 -> 5; 2 3, 3 2 -> 4
                end

                if m<=3 & n<=3
                    S(i,j,k,l)=Sr(m,n);
                elseif m>=3 & n>=3
                    S(i,j,k,l)=Sr(m,n)/4;
                else
                    S(i,j,k,l)=Sr(m,n)/2;
                end
            end
        end
    end
end