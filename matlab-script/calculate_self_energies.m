%% This file calculates the self-energy corrections for CNT bands

E_c=zeros(2,ikr_upper-ikr_lower+1);
E_v=zeros(2,ikr_upper-ikr_lower+1);
S_c=zeros(2,ikr_upper-ikr_lower+1);
S_v=zeros(2,ikr_upper-ikr_lower+1);
C_v=zeros(2,ikr_upper-ikr_lower+1,2);
C_c=zeros(2,ikr_upper-ikr_lower+1,2);

for ik=ikr_lower:ikr_upper
    ik
    mu_k=min_sub(1);
    k=mu_k*K1+ik*dk*K2;
    k_vec_tmp(ik-ikr_lower+1)=ik*dk;
    [Ek_v,Ek_c,Ck_v,Ck_c]=graphene_energy(e2p,t0,s0,a1,a2,k);
    E_v(1,ik-ikr_lower+1)=Ek_v;
    E_c(1,ik-ikr_lower+1)=Ek_c;
    C_v(1,ik-ikr_lower+1,:)=Ck_v;
    C_c(1,ik-ikr_lower+1,:)=Ck_c;

    for iq=-(nkc-1)/2:(nkc-1)/2
        for mu_q=1-Nu/2:Nu/2
            kp=k+mu_q*K1+iq*dk*K2;
            [Ekp_v,~,Ckp_v,Ckp_c]=graphene_energy(e2p,t0,s0,a1,a2,kp);
            vFT_d(:,:)=v_FT(mu_q+Nu,iq-iq_min+1,:,:);
            S_c(1,ik-ikr_lower+1)=S_c(1,ik-ikr_lower+1)- ...
                (conj(Ck_c(1))*Ckp_v(1)*vFT_d(1,1)*Ck_c(1)*conj(Ckp_v(1))+...
                conj(Ck_c(1))*Ckp_v(1)*vFT_d(1,2)*Ck_c(2)*conj(Ckp_v(2))+...
                conj(Ck_c(2))*Ckp_v(2)*vFT_d(2,1)*Ck_c(1)*conj(Ckp_v(1))+...
                conj(Ck_c(2))*Ckp_v(2)*vFT_d(2,2)*Ck_c(2)*conj(Ckp_v(2)))/(kappa*eps_q(mu_q+Nu,iq-iq_min+1));
            
            S_v(1,ik-ikr_lower+1)=S_v(1,ik-ikr_lower+1)- ...
                (conj(Ck_v(1))*Ckp_v(1)*vFT_d(1,1)*Ck_v(1)*conj(Ckp_v(1))+...
                conj(Ck_v(1))*Ckp_v(1)*vFT_d(1,2)*Ck_v(2)*conj(Ckp_v(2))+...
                conj(Ck_v(2))*Ckp_v(2)*vFT_d(2,1)*Ck_v(1)*conj(Ckp_v(1))+...
                conj(Ck_v(2))*Ckp_v(2)*vFT_d(2,2)*Ck_v(2)*conj(Ckp_v(2)))/(kappa*eps_q(mu_q+Nu,iq-iq_min+1));
        end;
    end;
    
    
    mu_k=min_sub(2);
    k=mu_k*K1+ik*dk*K2;
    
    [Ek_v,Ek_c,Ck_v,Ck_c]=graphene_energy(e2p,t0,s0,a1,a2,k);
    E_v(2,ik-ikr_lower+1)=Ek_v;
    E_c(2,ik-ikr_lower+1)=Ek_c;
    C_v(2,ik-ikr_lower+1,:)=Ck_v;
    C_c(2,ik-ikr_lower+1,:)=Ck_c;
    
    for iq=-(nkc-1)/2:(nkc-1)/2
        for mu_q=1-Nu/2:Nu/2
            kp=k+mu_q*K1+iq*dk*K2;
            [Ekp_v,~,Ckp_v,Ckp_c]=graphene_energy(e2p,t0,s0,a1,a2,kp);
            vFT_d(:,:)=v_FT(mu_q+Nu,iq-iq_min+1,:,:);
            
            S_c(2,ik-ikr_lower+1)=S_c(2,ik-ikr_lower+1)- ...
                (conj(Ck_c(1))*Ckp_v(1)*vFT_d(1,1)*Ck_c(1)*conj(Ckp_v(1))+...
                conj(Ck_c(1))*Ckp_v(1)*vFT_d(1,2)*Ck_c(2)*conj(Ckp_v(2))+...
                conj(Ck_c(2))*Ckp_v(2)*vFT_d(2,1)*Ck_c(1)*conj(Ckp_v(1))+...
                conj(Ck_c(2))*Ckp_v(2)*vFT_d(2,2)*Ck_c(2)*conj(Ckp_v(2)))/(kappa*eps_q(mu_q+Nu,iq-iq_min+1));
            
            S_v(2,ik-ikr_lower+1)=S_v(2,ik-ikr_lower+1)- ...
                (conj(Ck_v(1))*Ckp_v(1)*vFT_d(1,1)*Ck_v(1)*conj(Ckp_v(1))+...
                conj(Ck_v(1))*Ckp_v(1)*vFT_d(1,2)*Ck_v(2)*conj(Ckp_v(2))+...
                conj(Ck_v(2))*Ckp_v(2)*vFT_d(2,1)*Ck_v(1)*conj(Ckp_v(1))+...
                conj(Ck_v(2))*Ckp_v(2)*vFT_d(2,2)*Ck_v(2)*conj(Ckp_v(2)))/(kappa*eps_q(mu_q+Nu,iq-iq_min+1));
        end;
    end;
end;

%%
fig=fig+1; figure(fig); hold on;
plot(k_vec_tmp,E_v(1,:)/eV,'-r','LineWidth',2);
plot(k_vec_tmp,(E_c(1,:))/eV,'-b','LineWidth',2);
plot(k_vec_tmp,E_v(2,:)/eV,'--r','LineWidth',2);
plot(k_vec_tmp,(E_c(2,:))/eV,'--b','LineWidth',2);

plot(k_vec_tmp,(E_v(1,:)+real(S_v(1,:)))/eV,'-g','LineWidth',3);
plot(k_vec_tmp,(E_c(1,:)+real(S_c(1,:)))/eV,'-k','LineWidth',3);
plot(k_vec_tmp,(E_v(2,:)+real(S_v(2,:)))/eV,'--g','LineWidth',3);
plot(k_vec_tmp,(E_c(2,:)+real(S_c(2,:)))/eV,'--k','LineWidth',3);
axis tight; box on;

return;
%%
fig=fig+1; figure(fig); hold on;
plot(real(S_v(1,:))/eV,'-r','LineWidth',3);
plot(real(S_c(1,:))/eV,'-b','LineWidth',3);
plot(real(S_v(2,:))/eV,'--r','LineWidth',3);
plot(real(S_c(2,:))/eV,'--b','LineWidth',3);
axis tight; box on;