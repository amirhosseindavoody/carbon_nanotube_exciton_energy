%% This file calculates the direct and exchange term in the interaction kernel
%  based on the equation (15) in PRB 75,035407 (2007)
%  This file includes only the lowest subbands subbands of CNT as the basis.

ikr_min=-max(abs(ik_min-iKcm),abs(ik_max-iKcm)); %this is the lower limit of ikr (wave vector in reduced mass framework
ikr_max=+max(abs(ik_min-iKcm),abs(ik_max-iKcm)); %this is the upper limit of ikr (wave vector in reduced mass framework
nkr=ikr_max-ikr_min+1;

Kd=zeros(nkr*n_sub,nkr*n_sub);
Kx=zeros(nkr*n_sub,nkr*n_sub);
Ke=zeros(nkr*n_sub,nkr*n_sub); %this is the diagonal matrix containing quasiparticle transition energies


vFT_x(:,:)=v_FT(imin_sub,imin_sub,2*iKcm-iq_min+1,:,:); %"iq_min" is defined in calculate_dielectric_function.m

%% calculate intrasubband matrix elements between subband 1 and 1
imin_sub=1;
mu_kr=min_sub(imin_sub);
ipmin_sub=1;
mu_kpr=min_sub(ipmin_sub);
for ikr=ikr_min:ikr_max
    kr=mu_kr*K1+ikr*dk*K2;
    k_c=kr+Kcm;
    k_v=kr-Kcm;
    [~,Ek_c,~,Ck_c]=graphene_energy(e2p,t0,s0,a1,a2,k_c);
    [Ek_v,~,Ck_v,~]=graphene_energy(e2p,t0,s0,a1,a2,k_v);
    Ke((imin_sub-1)*nkr+(ikr-ikr_min+1),(imin_sub-1)*nkr+(ikr-ikr_min+1))=Ek_c-Ek_v;
    for ikpr=ikr_min:ikr_max
        kpr=mu_kpr*K1+ikpr*dk*K2;
        kp_c=kpr+Kcm;
        kp_v=kpr-Kcm;
        vFT_d(:,:)=v_FT(imin_sub,ipmin_sub,ikpr-ikr-iq_min+1,:,:); %"iq_min" is defined in calculate_dielectric_function.m
        
        [~,Ekp_c,~,Ckp_c]=graphene_energy(e2p,t0,s0,a1,a2,kp_c);
        [Ekp_v,~,Ckp_v,~]=graphene_energy(e2p,t0,s0,a1,a2,kp_v);
        
        Kd((imin_sub-1)*nkr+(ikr-ikr_min+1),(ipmin_sub-1)*nkr+(ikpr-ikr_min+1))= ...
            (conj(Ck_c(1))*Ckp_c(1)*vFT_d(1,1)*Ck_v(1)*conj(Ckp_v(1))+...
            conj(Ck_c(1))*Ckp_c(1)*vFT_d(1,2)*Ck_v(2)*conj(Ckp_v(2))+...
            conj(Ck_c(2))*Ckp_c(2)*vFT_d(2,1)*Ck_v(1)*conj(Ckp_v(1))+...
            conj(Ck_c(2))*Ckp_c(2)*vFT_d(2,2)*Ck_v(2)*conj(Ckp_v(2)))/(kappa*eps_q(imin_sub,ipmin_sub,ikpr-ikr-iq_min+1));
        
        Kx((imin_sub-1)*nkr+(ikr-ikr_min+1),(ipmin_sub-1)*nkr+(ikpr-ikr_min+1))= ...
            (conj(Ck_c(1))*Ck_v(1)*vFT_x(1,1)*Ckp_c(1)*conj(Ckp_v(1))+...
            conj(Ck_c(1))*Ck_v(1)*vFT_x(1,2)*Ckp_c(2)*conj(Ckp_v(2))+...
            conj(Ck_c(2))*Ck_v(2)*vFT_x(2,1)*Ckp_c(1)*conj(Ckp_v(1))+...
            conj(Ck_c(2))*Ck_v(2)*vFT_x(2,2)*Ckp_c(2)*conj(Ckp_v(2)));
    end;
end;

%% calculate intrasubband matrix elements between subband 2 and 2
imin_sub=2;
mu_kr=min_sub(imin_sub);
ipmin_sub=2;
mu_kpr=min_sub(ipmin_sub);
for ikr=ikr_max:-1:ikr_min
    kr=mu_kr*K1+ikr*dk*K2;
    k_c=kr+Kcm;
    k_v=kr-Kcm;
    [~,Ek_c,~,Ck_c]=graphene_energy(e2p,t0,s0,a1,a2,k_c);
    [Ek_v,~,Ck_v,~]=graphene_energy(e2p,t0,s0,a1,a2,k_v);
    Ke((imin_sub-1)*nkr+(ikr_max-ikr+1),(imin_sub-1)*nkr+(ikr_max-ikr+1))=Ek_c-Ek_v;
    for ikpr=ikr_max:-1:ikr_min
        kpr=mu_kpr*K1+ikpr*dk*K2;
        kp_c=kpr+Kcm;
        kp_v=kpr-Kcm;
        vFT_d(:,:)=v_FT(imin_sub,ipmin_sub,ikpr-ikr-iq_min+1,:,:); %"iq_min" is defined in calculate_dielectric_function.m
        
        [~,Ekp_c,~,Ckp_c]=graphene_energy(e2p,t0,s0,a1,a2,kp_c);
        [Ekp_v,~,Ckp_v,~]=graphene_energy(e2p,t0,s0,a1,a2,kp_v);
        
        Kd((imin_sub-1)*nkr+(ikr_max-ikr+1),(ipmin_sub-1)*nkr+(ikr_max-ikpr+1))= ...
            (conj(Ck_c(1))*Ckp_c(1)*vFT_d(1,1)*Ck_v(1)*conj(Ckp_v(1))+...
            conj(Ck_c(1))*Ckp_c(1)*vFT_d(1,2)*Ck_v(2)*conj(Ckp_v(2))+...
            conj(Ck_c(2))*Ckp_c(2)*vFT_d(2,1)*Ck_v(1)*conj(Ckp_v(1))+...
            conj(Ck_c(2))*Ckp_c(2)*vFT_d(2,2)*Ck_v(2)*conj(Ckp_v(2)))/(kappa*eps_q(imin_sub,ipmin_sub,ikpr-ikr-iq_min+1));
        
        Kx((imin_sub-1)*nkr+(ikr_max-ikr+1),(ipmin_sub-1)*nkr+(ikr_max-ikpr+1))= ...
            (conj(Ck_c(1))*Ck_v(1)*vFT_x(1,1)*Ckp_c(1)*conj(Ckp_v(1))+...
            conj(Ck_c(1))*Ck_v(1)*vFT_x(1,2)*Ckp_c(2)*conj(Ckp_v(2))+...
            conj(Ck_c(2))*Ck_v(2)*vFT_x(2,1)*Ckp_c(1)*conj(Ckp_v(1))+...
            conj(Ck_c(2))*Ck_v(2)*vFT_x(2,2)*Ckp_c(2)*conj(Ckp_v(2)));
    end;
end;

%% calculate intersubband matrix elements between subbands 1 and 2
imin_sub=1;
mu_kr=min_sub(imin_sub);
ipmin_sub=2;
mu_kpr=min_sub(ipmin_sub);
for ikr=ikr_min:ikr_max
    kr=mu_kr*K1+ikr*dk*K2;
    k_c=kr+Kcm;
    k_v=kr-Kcm;
    [~,Ek_c,~,Ck_c]=graphene_energy(e2p,t0,s0,a1,a2,k_c);
    [Ek_v,~,Ck_v,~]=graphene_energy(e2p,t0,s0,a1,a2,k_v);
    for ikpr=ikr_max:-1:ikr_min
        
        kpr=mu_kpr*K1+ikpr*dk*K2;
        kp_c=kpr+Kcm;
        kp_v=kpr-Kcm;
        vFT_d(:,:)=v_FT(imin_sub,ipmin_sub,ikpr-ikr-iq_min+1,:,:); %"iq_min" is defined in calculate_dielectric_function.m
        
        [~,Ekp_c,~,Ckp_c]=graphene_energy(e2p,t0,s0,a1,a2,kp_c);
        [Ekp_v,~,Ckp_v,~]=graphene_energy(e2p,t0,s0,a1,a2,kp_v);
        
        Kd((imin_sub-1)*nkr+(ikr-ikr_min+1),(ipmin_sub-1)*nkr+(ikr_max-ikpr+1))= ...
            (conj(Ck_c(1))*Ckp_c(1)*vFT_d(1,1)*Ck_v(1)*conj(Ckp_v(1))+...
            conj(Ck_c(1))*Ckp_c(1)*vFT_d(1,2)*Ck_v(2)*conj(Ckp_v(2))+...
            conj(Ck_c(2))*Ckp_c(2)*vFT_d(2,1)*Ck_v(1)*conj(Ckp_v(1))+...
            conj(Ck_c(2))*Ckp_c(2)*vFT_d(2,2)*Ck_v(2)*conj(Ckp_v(2)))/(kappa*eps_q(imin_sub,ipmin_sub,ikpr-ikr-iq_min+1));
        
        Kx((imin_sub-1)*nkr+(ikr-ikr_min+1),(ipmin_sub-1)*nkr+(ikr_max-ikpr+1))= ...
            (conj(Ck_c(1))*Ck_v(1)*vFT_x(1,1)*Ckp_c(1)*conj(Ckp_v(1))+...
            conj(Ck_c(1))*Ck_v(1)*vFT_x(1,2)*Ckp_c(2)*conj(Ckp_v(2))+...
            conj(Ck_c(2))*Ck_v(2)*vFT_x(2,1)*Ckp_c(1)*conj(Ckp_v(1))+...
            conj(Ck_c(2))*Ck_v(2)*vFT_x(2,2)*Ckp_c(2)*conj(Ckp_v(2)));
    end;
end;

%% calculate intersubband matrix elements between subbands 2 and 1
imin_sub=2;
mu_kr=min_sub(imin_sub);
ipmin_sub=1;
mu_kpr=min_sub(ipmin_sub);
for ikr=ikr_max:-1:ikr_min
    kr=mu_kr*K1+ikr*dk*K2;
    k_c=kr+Kcm;
    k_v=kr-Kcm;
    [~,Ek_c,~,Ck_c]=graphene_energy(e2p,t0,s0,a1,a2,k_c);
    [Ek_v,~,Ck_v,~]=graphene_energy(e2p,t0,s0,a1,a2,k_v);
    for ikpr=ikr_min:ikr_max
        
        kpr=mu_kpr*K1+ikpr*dk*K2;
        kp_c=kpr+Kcm;
        kp_v=kpr-Kcm;
        vFT_d(:,:)=v_FT(imin_sub,ipmin_sub,ikpr-ikr-iq_min+1,:,:); %"iq_min" is defined in calculate_dielectric_function.m
        
        [~,Ekp_c,~,Ckp_c]=graphene_energy(e2p,t0,s0,a1,a2,kp_c);
        [Ekp_v,~,Ckp_v,~]=graphene_energy(e2p,t0,s0,a1,a2,kp_v);
        
        Kd((imin_sub-1)*nkr+(ikr_max-ikr+1),(ipmin_sub-1)*nkr+(ikpr-ikr_min+1))= ...
            (conj(Ck_c(1))*Ckp_c(1)*vFT_d(1,1)*Ck_v(1)*conj(Ckp_v(1))+...
            conj(Ck_c(1))*Ckp_c(1)*vFT_d(1,2)*Ck_v(2)*conj(Ckp_v(2))+...
            conj(Ck_c(2))*Ckp_c(2)*vFT_d(2,1)*Ck_v(1)*conj(Ckp_v(1))+...
            conj(Ck_c(2))*Ckp_c(2)*vFT_d(2,2)*Ck_v(2)*conj(Ckp_v(2)))/(kappa*eps_q(imin_sub,ipmin_sub,ikpr-ikr-iq_min+1));
        
        Kx((imin_sub-1)*nkr+(ikr_max-ikr+1),(ipmin_sub-1)*nkr+(ikpr-ikr_min+1))= ...
            (conj(Ck_c(1))*Ck_v(1)*vFT_x(1,1)*Ckp_c(1)*conj(Ckp_v(1))+...
            conj(Ck_c(1))*Ck_v(1)*vFT_x(1,2)*Ckp_c(2)*conj(Ckp_v(2))+...
            conj(Ck_c(2))*Ck_v(2)*vFT_x(2,1)*Ckp_c(1)*conj(Ckp_v(1))+...
            conj(Ck_c(2))*Ck_v(2)*vFT_x(2,2)*Ckp_c(2)*conj(Ckp_v(2)));
    end;
end;

%% calculate free-particle transitions and make diagonal elements of the kernels real numbers.
Kd11=Kd(1:nkr,1:nkr);
Kd22=Kd(nkr+1:2*nkr,nkr+1:2*nkr);
Kd12=Kd(1:nkr,nkr+1:2*nkr);
Kd21=Kd(nkr+1:2*nkr,1:nkr);

Kx11=Kx(1:nkr,1:nkr);
Kx22=Kx(nkr+1:2*nkr,nkr+1:2*nkr);
Kx12=Kx(1:nkr,nkr+1:2*nkr);
Kx21=Kx(nkr+1:2*nkr,1:nkr);

Ke11=Ke(1:nkr,1:nkr);

Ef_tmp=diag(Ke11);

Kernel_intra=Ke11+2*Kx11-Kd11;
Kernel_inter=2*Kx12-Kd12;
Ex0_A1_tmp=eig(Kernel_intra-Kernel_inter);
Ex0_A2_tmp=eig(Kernel_intra+Kernel_inter);

Kernel_intra=Ke11-Kd11;
Kernel_inter=-Kd12;
Ex1_A1_tmp=eig(Kernel_intra-Kernel_inter);
Ex1_A2_tmp=eig(Kernel_intra+Kernel_inter);