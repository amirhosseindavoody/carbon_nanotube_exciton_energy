%% This file calculates the direct and exchange term in the interaction kernel
%  based on the equation (15) in PRB 75,035407 (2007)
%  This file includes only the lowest subbands subbands of CNT as the basis.

ikr_min=-max(abs(ik_min-iKcm),abs(ik_max-iKcm)); %this is the lower limit of ikr (wave vector in reduced mass framework
ikr_max=+max(abs(ik_min-iKcm),abs(ik_max-iKcm)); %this is the upper limit of ikr (wave vector in reduced mass framework
nkr=ikr_max-ikr_min+1;

Kd11=zeros(nkr,nkr);
Kd12=zeros(nkr,nkr);
Kx11=zeros(nkr,nkr);
Ke=zeros(nkr,nkr); %this is the diagonal matrix containing quasiparticle transition energies

Ck_c=zeros(2,1);
Ck_v=zeros(2,1);
Ckp_c=zeros(2,1);
Ckp_v=zeros(2,1);

vFT_x(:,:)=v_FT(Nu,2*iKcm-iq_min+1,:,:);

%% calculate intrasubband matrix elements between subband 1 and 1
imin_sub=1;
mu_kr=min_sub(imin_sub);
for ikr=ikr_min:ikr_max
    ikc=ikr+iKcm;
    Ek_c=E_c(1,ikc-ikr_lower+1);
    Sk_c=S_c(1,ikc-ikr_lower+1);
    Ck_c(:)=C_c(1,ikc-ikr_lower+1,:);
    ikv=ikr-iKcm;
    Ek_v=E_v(1,ikv-ikr_lower+1);
    Sk_v=S_v(1,ikv-ikr_lower+1);
    Ck_v(:)=C_v(1,ikv-ikr_lower+1,:);
    
    Ke((ikr-ikr_min+1),(ikr-ikr_min+1))=(Ek_c-Ek_v);
    Ke((ikr-ikr_min+1),(ikr-ikr_min+1))=Ke((ikr-ikr_min+1),(ikr-ikr_min+1))+real(Sk_c-Sk_v);
    
    %calculate intrasubband matrix elements between subbands 1 and 1
    ipmin_sub=1;
    mu_kpr=min_sub(ipmin_sub);
    for ikpr=ikr:ikr_max
        ikpc=ikpr+iKcm;
        Ckp_c(:)=C_c(1,ikpc-ikr_lower+1,:);
        ikpv=ikpr-iKcm;
        Ckp_v(:)=C_v(1,ikpv-ikr_lower+1,:);
        
        vFT_d(:,:)=v_FT(mu_kr-mu_kpr+Nu,ikr-ikpr-iq_min+1,:,:);
        Kd11((ikr-ikr_min+1),(ikpr-ikr_min+1))= ...
            (conj(Ck_c(1))*Ckp_c(1)*vFT_d(1,1)*Ck_v(1)*conj(Ckp_v(1))+...
            conj(Ck_c(1))*Ckp_c(1)*vFT_d(1,2)*Ck_v(2)*conj(Ckp_v(2))+...
            conj(Ck_c(2))*Ckp_c(2)*vFT_d(2,1)*Ck_v(1)*conj(Ckp_v(1))+...
            conj(Ck_c(2))*Ckp_c(2)*vFT_d(2,2)*Ck_v(2)*conj(Ckp_v(2)))/(kappa*eps_q(mu_kr-mu_kpr+Nu,ikr-ikpr-iq_min+1));

        Kx11((ikr-ikr_min+1),(ikpr-ikr_min+1))= ...
            (conj(Ck_c(1))*Ck_v(1)*vFT_x(1,1)*Ckp_c(1)*conj(Ckp_v(1))+...
            conj(Ck_c(1))*Ck_v(1)*vFT_x(1,2)*Ckp_c(2)*conj(Ckp_v(2))+...
            conj(Ck_c(2))*Ck_v(2)*vFT_x(2,1)*Ckp_c(1)*conj(Ckp_v(1))+...
            conj(Ck_c(2))*Ck_v(2)*vFT_x(2,2)*Ckp_c(2)*conj(Ckp_v(2)));
    end;
    
    %calculate intersubband matrix elements between subbands 1 and 2
    ipmin_sub=2;
    mu_kpr=min_sub(ipmin_sub);
    for ikpr=ikr_max:-1:-ikr
        ikpc=ikpr+iKcm;
        Ckp_c(:)=C_c(2,ikpc-ikr_lower+1,:);
        ikpv=ikpr-iKcm;
        Ckp_v(:)=C_v(2,ikpv-ikr_lower+1,:);
        
        vFT_d(:,:)=v_FT(mu_kr-mu_kpr+Nu,ikr-ikpr-iq_min+1,:,:);
        Kd12((ikr-ikr_min+1),(ikr_max-ikpr+1))= ...
            (conj(Ck_c(1))*Ckp_c(1)*vFT_d(1,1)*Ck_v(1)*conj(Ckp_v(1))+...
            conj(Ck_c(1))*Ckp_c(1)*vFT_d(1,2)*Ck_v(2)*conj(Ckp_v(2))+...
            conj(Ck_c(2))*Ckp_c(2)*vFT_d(2,1)*Ck_v(1)*conj(Ckp_v(1))+...
            conj(Ck_c(2))*Ckp_c(2)*vFT_d(2,2)*Ck_v(2)*conj(Ckp_v(2)))/(kappa*eps_q(mu_kr-mu_kpr+Nu,ikr-ikpr-iq_min+1));
    end;
end;

%% calculate free-particle transitions and make diagonal elements of the kernels real numbers.
Kd11=Kd11+ctranspose(Kd11);
Kx11=Kx11+ctranspose(Kx11);
Kd12=Kd12+ctranspose(Kd12);

Kd11=Kd11-diag(diag(Kd11))/2;
Kx11=Kx11-diag(diag(Kx11))/2;
Kd12=Kd12-diag(diag(Kd12))/2;

Ef_tmp=sort(diag(Ke));

K_A1=Ke-Kd11+Kd12;
K0_A2=Ke+4*Kx11-Kd11-Kd12;
K1_A2=Ke-Kd11-Kd12;

[Psi_A1_tmp,tmp]=eig(K_A1);
Ex_A1_tmp=diag(tmp);

[Psi0_A2_tmp,tmp]=eig(K0_A2);
Ex0_A2_tmp=diag(tmp);

[Psi1_A2_tmp,tmp]=eig(K1_A2);
Ex1_A2_tmp=diag(tmp);