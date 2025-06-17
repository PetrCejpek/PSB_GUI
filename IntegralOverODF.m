function Cav=IntegralOverODF(C,ODFinterpolant,sys)
% C - (full) tensor to average (3x3x3x3)
% ODFinterpolant - scattered interpolant of ODF (if ==1, we have an ideal
%                  polycrystal)
% sys - crystallographic system

[PHI1,PSI,PHI2]=meshgrid(0:5:360,0:5:180,0:5:360); % let's play with the grid in the future
dV=5*5*5*(pi/180)^3;

% ODF calculation
disp('Evaluating ODF...')
tic
switch sys
    case 'cubic'
         odfval=ODFinterpolant(PHI1,90-abs(PSI-90),mod(PHI2,90));
         %         odfval(m,:)=ODFinterpolant(Phi1(m,:),Psi(m,:),Phi2(m,:));
    case 'hexagonal'
         % IT NEEDS TO BE CHECKED
         odfval=ODFinterpolant(PHI1,90-abs(PSI-90),mod(PHI2,120));
end
toc

disp('Computing of rotated tensors...')
tic
Crot=zeros(3,3,3,3,size(PHI1,1),size(PHI1,2),size(PHI1,3)); % alokovat tak velké pole je problém, budu to muset integrovat krok po kroku
for n1=1:size(PHI1,1)
    for n2=1:size(PHI1,2)
        for n3=1:size(PHI1,3)
            R=RotMat(PHI2(n1,n2,n3),PSI(n1,n2,n3),PHI1(n1,n2,n3),'zxz'); % intrinsic rotation R=R1*R2*R3

            Crot(:,:,:,:,n1,n2,n3)=RotateFullElasticTensor(C,R);
        end
    end
end
toc

disp('Own integration...')
tic
Cav=zeros(size(C));
% Trapezoid integration in 3D
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
%                 disp([i,j,k,l])
                for n1=1:size(PHI1,1)-1
                    for n2=1:size(PHI1,2)-1
                        for n3=1:size(PHI1,3)-1
                            F_ef=(Crot(i,j,k,l,n1,  n2,  n3  ).*sind(PSI(n1,  n2,  n3  )).*odfval(n1,  n2,  n3)+...
                                  Crot(i,j,k,l,n1+1,n2,  n3  ).*sind(PSI(n1+1,n2,  n3  )).*odfval(n1+1,n2,  n3)+...
                                  Crot(i,j,k,l,n1+1,n2+1,n3  ).*sind(PSI(n1+1,n2+1,n3  )).*odfval(n1+1,n2+1,n3)+...
                                  Crot(i,j,k,l,n1,  n2+1,n3  ).*sind(PSI(n1,  n2+1,n3  )).*odfval(n1,  n2+1,n3)+...
                                  Crot(i,j,k,l,n1,  n2,  n3+1).*sind(PSI(n1,  n2,  n3+1)).*odfval(n1,  n2,  n3+1)+...
                                  Crot(i,j,k,l,n1+1,n2,  n3+1).*sind(PSI(n1+1,n2,  n3+1)).*odfval(n1+1,n2,  n3+1)+...
                                  Crot(i,j,k,l,n1+1,n2+1,n3+1).*sind(PSI(n1+1,n2+1,n3+1)).*odfval(n1+1,n2+1,n3+1)+...
                                  Crot(i,j,k,l,n1,  n2+1,n3+1).*sind(PSI(n1,  n2+1,n3+1)).*odfval(n1,  n2+1,n3+1))/8;

                            Cav(i,j,k,l)=Cav(i,j,k,l)+F_ef*dV;
                        end
                    end
                end
            end
        end
    end
end
toc

% normalisation term
disp('Computing normalisation...')
tic
A=0; % in fact, this integral should be equal to 8*pi^2 in the end, but depending on the choosen grid it may differ little bit
% Trapezoid integration in 3D
for n1=1:size(PHI1,1)-1
    for n2=1:size(PHI1,2)-1
        for n3=1:size(PHI1,3)-1
            F_ef=(sind(PSI(n1,  n2,  n3  )).*odfval(n1,  n2,  n3)+...
                  sind(PSI(n1+1,n2,  n3  )).*odfval(n1+1,n2,  n3)+...
                  sind(PSI(n1+1,n2+1,n3  )).*odfval(n1+1,n2+1,n3)+...
                  sind(PSI(n1,  n2+1,n3  )).*odfval(n1,  n2+1,n3)+...
                  sind(PSI(n1,  n2,  n3+1)).*odfval(n1,  n2,  n3+1)+...
                  sind(PSI(n1+1,n2,  n3+1)).*odfval(n1+1,n2,  n3+1)+...
                  sind(PSI(n1+1,n2+1,n3+1)).*odfval(n1+1,n2+1,n3+1)+...
                  sind(PSI(n1,  n2+1,n3+1)).*odfval(n1,  n2+1,n3+1))/8;
             A=A+F_ef*dV;
        end
    end
end
toc

Cav=Cav/A;