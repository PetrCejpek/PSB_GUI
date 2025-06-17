phi=0:1:359;
list1=cell(numel(phi),1);
list2=cell(numel(phi),1);
for k=1:numel(phi)
%     R=RotMat(0,0,0,'xyz'); % 100
%     R=RotMat(45,06,0,'xyz'); % 110
    R=RotMat(45,54.736,0,'zxz'); % 111
    Rphi=RotMat(0,0,phi(k),'xyz');
    list{k}=Rphi*R;
    vec1=Rphi*R*[-0.5;-0.5;-0.5];
    vec2=Rphi*R*[0.5;-0.5;-0.5];
    vec3=Rphi*R*[0.5;0.5;-0.5];
    vec4=Rphi*R*[-0.5;0.5;-0.5];
    vec5=Rphi*R*[-0.5;-0.5;0.5];
    vec6=Rphi*R*[0.5;-0.5;0.5];
    vec7=Rphi*R*[0.5;0.5;0.5];
    vec8=Rphi*R*[-0.5;0.5;0.5];
    list2{k}=[vec1 vec2 vec3 vec4 vec5 vec6 vec7 vec8]';
end

figure
hold on
for k=1:numel(phi)
%     list{k}=Rphi*R;
%     list2{k}=[vec1 vec2 vec3 vec4 vec5 vec6 vec7 vec8]';

    plot3([list2{k}(1,1) list2{k}(2,1) list2{k}(3,1) list2{k}(4,1) list2{k}(1,1)],...
        [list2{k}(1,2) list2{k}(2,2) list2{k}(3,2) list2{k}(4,2) list2{k}(1,2)],...
        [list2{k}(1,3) list2{k}(2,3) list2{k}(3,3) list2{k}(4,3) list2{k}(1,3)],'b')
    plot3([list2{k}(5,1) list2{k}(6,1) list2{k}(7,1) list2{k}(8,1) list2{k}(5,1)],...
        [list2{k}(5,2) list2{k}(6,2) list2{k}(7,2) list2{k}(8,2) list2{k}(5,2)],...
        [list2{k}(5,3) list2{k}(6,3) list2{k}(7,3) list2{k}(8,3) list2{k}(5,3)],'b')
    plot3([list2{k}(1,1) list2{k}(5,1)],[list2{k}(1,2) list2{k}(5,2)],[list2{k}(1,3) list2{k}(5,3)],'b')
    plot3([list2{k}(2,1) list2{k}(6,1)],[list2{k}(2,2) list2{k}(6,2)],[list2{k}(2,3) list2{k}(6,3)],'b')
    plot3([list2{k}(3,1) list2{k}(7,1)],[list2{k}(3,2) list2{k}(7,2)],[list2{k}(3,3) list2{k}(7,3)],'b')
    plot3([list2{k}(4,1) list2{k}(8,1)],[list2{k}(4,2) list2{k}(8,2)],[list2{k}(4,3) list2{k}(8,3)],'b')

%     % 110
%     plot3([0 list2{k}(8,1)],[0 list2{k}(8,2)],[0 list2{k}(8,3)],'r')
%     plot3([0 list2{k}(5,1)+list2{k}(6,1)]/2,[0 list2{k}(5,2)+list2{k}(6,2)]/2,[0 list2{k}(5,3)+list2{k}(6,3)]/2,'g')
    
    % 111
    plot3([0 list2{k}(6,1)],[0 list2{k}(6,2)],[0 list2{k}(6,3)],'r')
    plot3([0 list2{k}(5,1)+list2{k}(1,1)]/2,[0 list2{k}(5,2)+list2{k}(1,2)]/2,[0 list2{k}(5,3)+list2{k}(1,3)]/2,'g')
end

hold off
axis equal
xlabel('x')
ylabel('y')
zlabel('z')