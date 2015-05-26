%% draw graphene reciprocal lattice unit cell

tmp1 = (b1-b2)/3;
Rot=[cos(2*pi/3),sin(2*pi/3);-sin(2*pi/3),cos(2*pi/3)];
tmp2 = transpose(Rot*transpose(tmp1));
tmp3 = transpose(Rot*transpose(tmp2));

K1=(-t2*b1+t1*b2)/Nu;
K2=(mC*b1-nC*b2)/Nu;

fig=fig+1; figure(fig); hold on; box on;
quiver(0,0,b1(1),b1(2),'AutoScale','off','LineWidth',3);
quiver(0,0,b2(1),b2(2),'AutoScale','off','LineWidth',3);
% quiver(0,0,tmp1(1),tmp1(2),'AutoScale','off','LineWidth',3);
% quiver(0,0,tmp2(1),tmp2(2),'AutoScale','off','LineWidth',3);
% quiver(0,0,tmp3(1),tmp3(2),'AutoScale','off','LineWidth',3);

tmpScale=5;
quiver(0,0,tmpScale*K1(1),tmpScale*K1(2),'AutoScale','off','LineWidth',3);
quiver(0,0,tmpScale*K2(1),tmpScale*K2(2),'AutoScale','off','LineWidth',3);

tmpScale=10;
for i=-10:10
    plot([0+i*K1(1),tmpScale*K2(1)+i*K1(1)],[0+i*K1(2),tmpScale*K2(2)+i*K1(2)],'-k','LineWidth',2);
end;
axis equal; axis tight;

brill_zone_coor = zeros(7,2);
brill_zone_coor(1,:) = [0,0];
brill_zone_coor(2,:) = brill_zone_coor(1,:) + tmp1;
brill_zone_coor(3,:) = brill_zone_coor(2,:) - tmp3;
brill_zone_coor(4,:) = brill_zone_coor(3,:) + tmp2;
brill_zone_coor(5,:) = brill_zone_coor(4,:) - tmp1;
brill_zone_coor(6,:) = brill_zone_coor(5,:) + tmp3;
brill_zone_coor(7,:) = brill_zone_coor(6,:) - tmp2;

brill_zone_coor(:,1) = brill_zone_coor(:,1) + tmp3(1);
brill_zone_coor(:,2) = brill_zone_coor(:,2) + tmp3(2);

% fig=fig+1; figure(fig); hold on; box on;
plot(brill_zone_coor(:,1),brill_zone_coor(:,2),'LineWidth',3);