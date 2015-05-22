%% This file calculates the geometrical properties of CNT and the
% coordinates of atoms in the CNT unit cell.

%% Unit vectors and reciprocal lattice vectors
a1=[sqrt(3)/2,1/2]*aL;
a2=[sqrt(3)/2,-1/2]*aL;
b1=[1/sqrt(3),1]*2*pi/aL;
b2=[1/sqrt(3),-1]*2*pi/aL;

%% Calculate chirality and translational vectors
ch_vec=nC*a1+mC*a2; %chiral vector
len_ch=aL*sqrt(nC^2+mC^2+nC*mC); %length of the chiral vector
radius=len_ch/2/pi;

dR=gcd(2*nC+mC,2*mC+nC);
t1=(2*mC+nC)/dR;
t2=-(2*nC+mC)/dR;
t_vec=t1*a1+t2*a2; %translational vector
Nu=2*(nC^2+mC^2+nC*mC)/dR; %number of hexagons in CNT unit cell

p1=max([(Nu/nC+1/t1)/(mC/nC-t2/t1),(1/nC+1/t1)/(mC/nC-t2/t1)]);
p2=min([(Nu/nC+1/t1)/(mC/nC-t2/t1),(1/nC+1/t1)/(mC/nC-t2/t1)]);
for i=ceil(p2):floor(p1)
    p=i;
    q=t2/t1*p+1/t1;
    if (q==ceil(q))
        break;
    elseif(i==floor(p1))
        disp('MC not found!');
        exit;
    end;
end;
tmp=gcd(p,q);
p=p/tmp;
q=q/tmp;
MC=mC*p-nC*q;

aCC_vec=1/3*(a1+a2);

%% Calculate coordinate of lattice sites in the 2D unit cell.
k=0;
pos=zeros(1,2);
posA=zeros(1,2);
posB=zeros(1,2);
for i=0:t1+nC
    for j=t2:mC
        if(((t2/t1*i-j)<=0)&&((mC/nC*i-j)>=0) ...
                &&((t2/t1*(i-nC)-(j-mC))>0)&&((mC/nC*(i-t1)-(j-t2))<0))
            k=k+1;
            pos(k,:)=i*a1+j*a2;
            posA(k,:)=pos(k,:);
            posB(k,:)=pos(k,:)+aCC_vec;
        end;
    end;
end;

% here we rotate the coordinates of atoms, chiral vector, and translational
% vectors but we do not rotate lattice vectors (a1 and a1) and we do NOT
% rotate reciprocal lattice vectors (b1 and b2).
cosTh=ch_vec(1)/(norm(ch_vec));
sinTh=ch_vec(2)/(norm(ch_vec));
Rot=[cosTh,sinTh;-sinTh,cosTh];
ch_vec=transpose(Rot*transpose(ch_vec));
t_vec=transpose(Rot*transpose(t_vec));
a1=transpose(Rot*transpose(a1));
a2=transpose(Rot*transpose(a2));
b1=transpose(Rot*transpose(b1));
b2=transpose(Rot*transpose(b2));
aCC_vec=transpose(Rot*transpose(aCC_vec));

for k=1:Nu
    tmp(1,1)=pos(k,1);
    tmp(2,1)=pos(k,2);
    pos(k,:)=Rot*tmp;
    tmp(1,1)=posA(k,1);
    tmp(2,1)=posA(k,2);
    posA(k,:)=Rot*tmp;
    tmp(1,1)=posB(k,1);
    tmp(2,1)=posB(k,2);
    posB(k,:)=Rot*tmp;
    if (posB(k,1)>ch_vec(1))
        posB(k,1)=posB(k,1)-ch_vec(1);
    end;
    if (posB(k,2)<0)
        posB(k,2)=posB(k,2)+t_vec(2);
    end;
end;

% these are the distance between atoms of type A and B in a CNT unit cell
% considering the circular symmetry of the real CNT structure.
posAA=zeros(1,2);
posAB=zeros(1,2);
posBA=zeros(1,2);
posBB=zeros(1,2);
for k=1:Nu
    posAA(k,:)=posA(k,:)-posA(1,:);
    posAB(k,:)=posA(k,:)-posB(1,:);
    posBA(k,:)=posB(k,:)-posA(1,:);
    posBB(k,:)=posB(k,:)-posB(1,:);
    if (posAA(k,1)>(ch_vec(1)/2))
        posAA(k,1)=posAA(k,1)-ch_vec(1);
    end;
    if (posAB(k,1)>(ch_vec(1)/2))
        posAB(k,1)=posAB(k,1)-ch_vec(1);
    end;
    if (posBA(k,1)>(ch_vec(1)/2))
        posBA(k,1)=posBA(k,1)-ch_vec(1);
    end;
    if (posBB(k,1)>(ch_vec(1)/2))
        posBB(k,1)=posBB(k,1)-ch_vec(1);
    end;
end;

clear tmp;

posTOT=zeros(1,2*Nu); %this variable stores just the x-coordinate of all the atoms in unit cell for the calculation of tight-binding coefficients.
for k=1:Nu
    posTOT(1,k)=posA(k,1);
    posTOT(1,Nu+k)=posB(k,1);
end;

% %plot CNT unit cell
fig=fig+1; figure(fig); hold on; box on;
plot(posA(:,1),posA(:,2),'k.','MarkerSize',20);
plot(posB(:,1),posB(:,2),'g.','MarkerSize',20);

plot(posA(1,1),posA(1,2),'b.','MarkerSize',20);
plot(posB(1,1),posB(1,2),'r.','MarkerSize',20);

ch_plot=zeros(2,2);
ch_plot(2,:)=ch_vec;
plot(ch_plot(:,1),ch_plot(:,2),'-b','LineWidth',3);

t_plot=zeros(2,2);
t_plot(2,:)=t_vec;
plot(t_plot(:,1),t_plot(:,2),'-r','LineWidth',3);

ch_plot(1,:)=t_vec;
ch_plot(2,:)=ch_vec+t_vec;
plot(ch_plot(:,1),ch_plot(:,2),'--b','LineWidth',3);

t_plot(1,:)=ch_vec;
t_plot(2,:)=ch_vec+t_vec;
plot(t_plot(:,1),t_plot(:,2),'--r','LineWidth',3);
axis equal; axis tight;

fig=fig+1; figure(fig); hold on; box on;
a1_plot=zeros(2,2);
a1_plot(2,:)=a1/norm(a1);
plot(a1_plot(:,1),a1_plot(:,2),'-b','LineWidth',3);

a1_plot=zeros(2,2);
a1_plot(2,:)=a2/norm(a2);
plot(a1_plot(:,1),a1_plot(:,2),'-r','LineWidth',3);

a1_plot=zeros(2,2);
a1_plot(2,:)=b1/norm(b1);
plot(a1_plot(:,1),a1_plot(:,2),'--b','LineWidth',3);

a1_plot=zeros(2,2);
a1_plot(2,:)=b2/norm(b2);
plot(a1_plot(:,1),a1_plot(:,2),'--r','LineWidth',3);
axis equal; axis tight;

%% Calculate coordinates of atoms in 3D unit cell
A_coor=zeros(Nu,3);
tmp=2*pi*posA(:,1)/len_ch;
A_coor(:,1)=radius*sin(tmp);
A_coor(:,2)=posA(:,2);
A_coor(:,3)=-radius*cos(tmp);

B_coor=zeros(Nu,3);
tmp=2*pi*posB(:,1)/len_ch;
B_coor(:,1)=radius*sin(tmp);
B_coor(:,2)=posB(:,2);
B_coor(:,3)=-radius*cos(tmp);

clear tmp;

% fig=fig+1; figure(fig); hold on; box on;
% plot3(A_coor(:,1),A_coor(:,2),A_coor(:,3),'k.','MarkerSize',20);
% plot3(B_coor(:,1),B_coor(:,2),B_coor(:,3),'k.','MarkerSize',20);
% axis equal;