%% This file plots electron-hole amplitude which can be thought of as the exciton wavefunction

ikr_min=-max(abs(ik_min-iKcm),abs(ik_max-iKcm)); %this is the lower limit of ikr (wave vector in reduced mass framework
ikr_max=+max(abs(ik_min-iKcm),abs(ik_max-iKcm)); %this is the upper limit of ikr (wave vector in reduced mass framework
nkr=ikr_max-ikr_min+1;

Psi0_A2=Psi0_A2_tmp;
xstate=1;

deltaR=[0,0];
deltaC=[0,0,0];
k=0;
for u=-2:1
    for s1=1:Nu
        k=k+1;
        deltaR(k,:)=u*t_vec+posAA(s1,:);
        deltaC(k,:)=u*[t_vec,0]+A_coor(s1,:);
    end;
end

Psi_real=zeros(k,1);

mu_kr=min_sub(1);
for ikr=ikr_min:ikr_max
    kr=mu_kr*K1+ikr*dk*K2;
    k_c=kr+Kcm;
    k_v=kr-Kcm;
    [~,Ek_c1,~,Ck_c1]=graphene_energy(e2p,t0,s0,a1,a2,k_c);
    [Ek_v1,~,Ck_v1,~]=graphene_energy(e2p,t0,s0,a1,a2,k_v);
    
    k_c=-kr+Kcm;
    k_v=-kr-Kcm;
    [~,Ek_c2,~,Ck_c2]=graphene_energy(e2p,t0,s0,a1,a2,k_c);
    [Ek_v2,~,Ck_v2,~]=graphene_energy(e2p,t0,s0,a1,a2,k_v);
    
    Psi_real=Psi_real+Psi0_A2(ikr-ikr_min+1,xstate)*(Ck_c1(1)*conj(Ck_v1(1))*exp(1i*(kr(1)*deltaR(:,1)+kr(2)*deltaR(:,2)))+Ck_c2(1)*conj(Ck_v2(1))*exp(-1i*(kr(1)*deltaR(:,1)+kr(2)*deltaR(:,2))));
end;

fig=fig+1; figure(fig); hold on;
scatter(deltaR(:,2),deltaR(:,1),abs(1e1*Psi_real(:,1)),abs(1e1*Psi_real(:,1)),'Fill');
axis tight; box on;
pbaspect([10,1,1]);
% daspect ([1 1 1]);

fig=fig+1; figure(fig); hold on;
scatter3(deltaC(:,1),deltaC(:,2),deltaC(:,3),abs(2e1*Psi_real(:,1)),abs(Psi_real(:,1)),'Fill');
pbaspect('manual');
pbaspect([1,10,1]);
% daspect ([1 1 1]);
axis tight; box on;
rotate3d;
axis(gca,'vis3d');

return;
%%
figure
plot(imag(Psi_real),'LineWidth',3);
return;
%%
figure; hold on;
kr_vec=(ikr_min:ikr_max)*dk;
plot(kr_vec,real(Psi0_A2(:,2)));
plot(kr_vec,-real(Psi0_A2(:,1)));