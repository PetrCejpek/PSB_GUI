function y=Reci(a,b,c,al,be,ga)
% vraci matici s reciprokymi vektory
% a,b,c ... mrizove parametry
% al,be,ga ... uhly v elementarni bunce

% vektory prime mrize, vektor a miri podel osy x
% a_p=[a; 0; 0];
% b_p=[b*cos(ga*pi/180); b*sin(ga*pi/180); 0];
% c_p(1,1)=c*cos(be*pi/180);
% c_p(2,1)=c*(cos(al*pi/180)-cos(ga*pi/180)*cos(be*pi/180))/sin(ga*pi/180);
% c_p(3,1)=sqrt(c^2-c_p(1,1)^2-c_p(2,1)^2);
A=Prima(a,b,c,al,be,ga);

y=2*pi*inv(A)';