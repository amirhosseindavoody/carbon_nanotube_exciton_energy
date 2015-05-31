%% This file visualizes the results of the fortran program for CNT bethe salpeter equation
clear all; clc; fig=10;
close all;
dir='C:\Users\amirhossein\Google Drive\Research\Exciton\Data\Environmental Effect\CNT-Exciton\CNT(08,00)-nkg(1001)-nr(0200)-E_th(1.5)-Kcm_max(1.5)-i_sub(1)-Ckappa(2.0)\';
eV=1.6e-19;

%% plot CNT unit cell
FileName=[dir,'posA.dat'];
posA=load(FileName);
FileName=[dir,'posB.dat'];
posB=load(FileName);

fig=fig+1; figure(fig); hold on; box on;
plot(posA(:,1),posA(:,2),'.','MarkerSize',20);
plot(posB(:,1),posB(:,2),'.','MarkerSize',20);

% plot(posA(1,1),posA(1,2),'g.','MarkerSize',20);
% plot(posB(1,1),posB(1,2),'k.','MarkerSize',20);

axis equal; axis tight;

%% plot CNT energy dispersion
FileName=[dir,'CondBand.dat'];
Ec_tmp=load(FileName);
FileName=[dir,'ValeBand.dat'];
Ev_tmp=load(FileName);

[Nu,nkc]=size(Ec_tmp);
Nu=Nu-1;
k_vec=Ec_tmp(1,:);

E_c=Ec_tmp(2:Nu+1,:);
E_v=Ev_tmp(2:Nu+1,:);

fig=fig+1; figure(fig); hold on; box on;
plot(k_vec,E_c/eV,'-','LineWidth',3);
plot(k_vec,E_v/eV,'-','LineWidth',3);
axis tight;

%% plot energy dispersion of the considered subbands
FileName=[dir,'CondBand_Sub.dat'];
Ec_tmp=load(FileName);
FileName=[dir,'ValeBand_Sub.dat'];
Ev_tmp=load(FileName);

[tmp,nkr]=size(Ec_tmp);
tmp=tmp-1;
k_vec=Ec_tmp(1,:);

E_c=Ec_tmp(2:tmp+1,:);
E_v=Ev_tmp(2:tmp+1,:);

figure(fig); hold on; box on;
plot(k_vec,E_c/eV,'-','LineWidth',3);
plot(k_vec,E_v/eV,'-','LineWidth',3);
axis tight;

%% plot self-energy corrections
FileName=[dir,'CondSelfEnergy_Sub.dat'];
S_c=load(FileName);
FileName=[dir,'ValeSelfEnergy_Sub.dat'];
S_v=load(FileName);

fig=fig+1; figure(fig); hold on; box on;
plot(k_vec,E_c/eV,'-','LineWidth',2);
plot(k_vec,E_v/eV,'-','LineWidth',2);
plot(k_vec,(E_c+S_c)/eV,'-','LineWidth',3);
plot(k_vec,(E_v+S_v)/eV,'-','LineWidth',3);
axis tight;

%% plot eps(0,q)
dk=k_vec(2)-k_vec(1);
FileName=[dir,'eps_q.dat'];
eps_q=load(FileName);
[n_mu,n_q]=size(eps_q);
mu_vec= -Nu+1:Nu-1;

q_max=floor(n_q/2)*dk;
q_min=-floor(n_q/2)*dk;
q_vec=linspace(q_min,q_max,n_q);

fig=fig+1; figure(fig); hold on; box on;
plot(q_vec,(eps_q(Nu,:)),'-','LineWidth',3);
axis tight;
% return;

%% plot exciton energy Ex0_Ep
% FileName=[dir,'Ex0_Ep.dat'];
% Ex0_Ep=load(FileName);
% [nKcm,nX]=size(Ex0_Ep);
% Kcm_vec=dk*(-(nKcm-1)/2:+(nKcm-1)/2);
% 
% fig=fig+1; figure(fig); hold on; box on;
% for i=1:nX
%     plot(Kcm_vec,Ex0_Ep(:,i)/eV,'-k','LineWidth',3);
% end;
% axis tight;
% 
% %% plot exciton energy Ex0_Em
% FileName=[dir,'Ex0_Em.dat'];
% Ex0_Em=load(FileName);
% [nKcm,nX]=size(Ex0_Em);
% Kcm_vec=dk*(-(nKcm-1)/2:+(nKcm-1)/2);
% 
% % fig=fig+1; figure(fig); hold on; box on;
% for i=1:nX
%     plot(Kcm_vec,Ex0_Em(:,i)/eV,'-b','LineWidth',3);
% end;
% axis tight;
% 
% 
% %% plot exciton energy Ex1_Ep
% FileName=[dir,'Ex1_Ep.dat'];
% Ex1_Ep=load(FileName);
% [nKcm,nX]=size(Ex1_Ep);
% Kcm_vec=dk*(-(nKcm-1)/2:+(nKcm-1)/2);
% 
% % fig=fig+1; figure(fig); hold on; box on;
% for i=1:nX
%     plot(Kcm_vec,Ex1_Ep(:,i)/eV,'-k','LineWidth',3);
% end;
% axis tight;
% 
% %% plot exciton energy Ex1_Em
% FileName=[dir,'Ex1_Em.dat'];
% Ex1_Em=load(FileName);
% [nKcm,nX]=size(Ex1_Em);
% Kcm_vec=dk*(-(nKcm-1)/2:+(nKcm-1)/2);
% 
% % fig=fig+1; figure(fig); hold on; box on;
% for i=1:nX
%     plot(Kcm_vec,Ex1_Em(:,i)/eV,'-g','LineWidth',3);
% end;
% axis tight;

% return;

%% plot exciton energy Ex_A1
FileName=[dir,'Ex_A1.dat'];
Ex_A1=load(FileName);
[nKcm,nX]=size(Ex_A1);
Kcm_vec=dk*(-(nKcm-1)/2:+(nKcm-1)/2);

fig=fig+1; figure(fig); hold on; box on;
for i=1:nX
    plot(Kcm_vec,Ex_A1(:,i)/eV,'-','LineWidth',3);
end;
axis tight;

%% plot exciton energy Ex0_A2
FileName=[dir,'Ex0_A2.dat'];
Ex0_A2=load(FileName);
[nKcm,nX]=size(Ex0_A2);
Kcm_vec=dk*(-(nKcm-1)/2:+(nKcm-1)/2);

fig=fig+1; figure(fig); hold on; box on;
for i=1:nX
    plot(Kcm_vec,Ex0_A2(:,i)/eV,'-','LineWidth',3);
end;
axis tight;

%% plot exciton energy Ex_A1
FileName=[dir,'Ex1_A2.dat'];
Ex1_A2=load(FileName);
[nKcm,nX]=size(Ex1_A2);
Kcm_vec=dk*(-(nKcm-1)/2:+(nKcm-1)/2);

fig=fig+1; figure(fig); hold on; box on;
for i=1:nX
    plot(Kcm_vec,Ex1_A2(:,i)/eV,'-','LineWidth',3);
end;
axis tight;
return;

%% plot exciton wavefunction in k-space
FileName=[dir,'Psi_A1.dat'];
tmp=load(FileName);
nKcm=size(tmp,1);
nk=size(tmp,2)/nX/2;
Psi_A1=zeros(nKcm,nX,nk);

for i=1:nX
    for j=1:nk
        Psi_A1(:,i,j)=tmp(:,(i-1)*2*nk+2*j-1)+1i*tmp(:,(i-1)*2*nk+2*j);
    end;
end;
clear tmp;

FileName=[dir,'Psi0_A2.dat'];
tmp=load(FileName);
[nKcm,ncol]=size(tmp);
nk=ncol/nX/2;
Psi0_A2=zeros(nKcm,nX,nk);

for i=1:nX
    for j=1:nk
        Psi0_A2(:,i,j)=tmp(:,(i-1)*2*nk+2*j-1)+1i*tmp(:,(i-1)*2*nk+2*j);
    end;
end;
clear tmp;

FileName=[dir,'Psi1_A2.dat'];
tmp=load(FileName);
[nKcm,ncol]=size(tmp);
nk=ncol/nX/2;
kr_vec=dk*(-(nk-1)/2:+(nk-1)/2);
Psi1_A2=zeros(nKcm,nX,nk);

for i=1:nX
    for j=1:nk
        Psi1_A2(:,i,j)=tmp(:,(i-1)*2*nk+2*j-1)+1i*tmp(:,(i-1)*2*nk+2*j);
    end;
end;
clear tmp;

iKcm=(nKcm+1)/2;
fig=fig+1; figure(fig); box on; hold on;
for iX=1
%     tmp(1,:)=Psi_A1(iKcm,iX,:);
%     plot(kr_vec,abs(tmp),'-','LineWidth',3);
    tmp(1,:)=Psi0_A2(iKcm,iX,:);
    plot(kr_vec,abs(tmp),'-','LineWidth',3);
%     tmp(1,:)=Psi1_A2(iKcm,iX,:);
%     plot(kr_vec,abs(tmp),'-','LineWidth',3);
end;
axis tight;
clear tmp;