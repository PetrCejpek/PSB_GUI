function Trot=RotateFullElasticTensor(T,R)
% T is the 3x3x3x3 tensor/matrix
% Tr is the 6x6 tensor/matrix

Trot=zeros(3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                % sumation indices
                for m=1:3
                    for n=1:3
                        for o=1:3
                            for p=1:3
                                Trot(i,j,k,l)=Trot(i,j,k,l)+R(m,i)*R(n,j)*R(o,k)*R(p,l)*T(m,n,o,p);
                            end
                        end
                    end
                end
            end
        end
    end
end