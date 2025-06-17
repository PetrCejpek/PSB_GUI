function Tr=Tensor4_rotation(T,R)
% Function Tensor4_rotation will rotate 4th order tensor (size 3x3x3x3)
% with the use of a rotation matrix R (size 3x3)

Tr=zeros(size(T));
for i2=1:3
    for j2=1:3
        for k2=1:3
            for l2=1:3
                Tr(i2,j2,k2,l2)=0;

                for i1=1:3
                    for j1=1:3
                        for k1=1:3
                            for l1=1:3
                                Tr(i2,j2,k2,l2)=Tr(i2,j2,k2,l2)+Tr(i1,j1,k1,l1)*R(i2,i1)*R(j2,j1)*R(k2,k1)*R(l2,l1);
                            end
                        end
                    end
                end
            end
        end
    end
end