%% This file calculates dielectric conductivity and Fourier transform of
% Coulomb interaction based on the equations (13) and (17) in PRB 75,035407 (2007)

% It should be noted that coefficients C_c and C_v are not periodic over
% the first birillouin zone of CNT.

%% Calculate PI(q) from equation (17): PRB 75,035407 (2007)
PI=zeros(2*Nu-1,iq_max-iq_min+1);

for ik=-(nkc-1)/2:(nkc-1)/2
    ik
    for mu_k=1-Nu/2:Nu/2
        k=mu_k*K1+ik*dk*K2;
        [Ek_v,Ek_c,Ck_v,Ck_c]=graphene_energy(e2p,t0,s0,a1,a2,k);
        for iq=iq_min:iq_max
            for mu=(1-Nu):(Nu-1)
                mu_kq=mu_k+mu;
                kq=mu_kq*K1+(ik+iq)*dk*K2;
                [Ekq_v,Ekq_c,Ckq_v,Ckq_c]=graphene_energy(e2p,t0,s0,a1,a2,kq);
                PI(mu+Nu,iq-iq_min+1)=PI(mu+Nu,iq-iq_min+1)+ ...
                    (abs(ctranspose(Ck_v)*Ckq_c))^2/(Ekq_c-Ek_v)+(abs(ctranspose(Ck_c)*Ckq_v))^2/(Ek_c-Ekq_v);
            end;
        end;
    end;
end;
% return;
%% Calculate the Fourier transform of Coulomb potential from equation (13)

% the commented code below this section has a slightly higher value which
% matches the result in the reference paper better. you have to check the cause of
% this mismatch later.
v_FT=zeros(2*Nu-1,iq_max-iq_min+1,2,2);
v_q=zeros(2*Nu-1,iq_max-iq_min+1);
tmp=zeros(1,iq_max-iq_min+1,1,1);
tmp10=1;

for iq=iq_min:iq_max
    q_vec(iq-iq_min+1)=iq*dk;
end;

for mu=(1-Nu):(Nu-1)
    mu
    for u=-floor(2*nr):floor(2*nr)
        for s1=1:Nu
            deltaR=u*t_vec+posAA(s1,:);
            tmp(1,:,1,1)=(exp(1i*(q_vec*(K2*transpose(deltaR))+mu*(K1*transpose(deltaR))))*Upp/(sqrt((4*pi*eps0/q0^2*Upp*norm(deltaR))^2+1)))';
            v_FT(mu+Nu,:,1,1)=v_FT(mu+Nu,:,1,1)+tmp;
            
            deltaR=u*t_vec+posBA(s1,:);
            tmp(1,:,1,1)=(exp(1i*(q_vec*(K2*transpose(deltaR))+mu*(K1*transpose(deltaR))))*Upp/(sqrt((4*pi*eps0/q0^2*Upp*norm(deltaR))^2+1)))';
            v_FT(mu+Nu,:,1,2)=v_FT(mu+Nu,:,1,2)+tmp;
            
            deltaR=u*t_vec+posAB(s1,:);
            tmp(1,:,1,1)=(exp(1i*(q_vec*(K2*transpose(deltaR))+mu*(K1*transpose(deltaR))))*Upp/(sqrt((4*pi*eps0/q0^2*Upp*norm(deltaR))^2+1)))';
            v_FT(mu+Nu,:,2,1)=v_FT(mu+Nu,:,2,1)+tmp;
            
            deltaR=u*t_vec+posBB(s1,:);
            tmp(1,:,1,1)=(exp(1i*(q_vec*(K2*transpose(deltaR))+mu*(K1*transpose(deltaR))))*Upp/(sqrt((4*pi*eps0/q0^2*Upp*norm(deltaR))^2+1)))';
            v_FT(mu+Nu,:,2,2)=v_FT(mu+Nu,:,2,2)+tmp;
        end;
    end;
    for i=iq_min:iq_max
        v_q(mu+Nu,i-iq_min+1)=1/4*sum(sum(v_FT(mu+Nu,i-iq_min+1,:,:)));
    end;
end;

v_FT=v_FT/(nkc*Nu);
v_q=v_q/(nkc*Nu);

eps_q=1+v_q.*PI;

%% Calculate the Fourier transform of Coulomb potential from equation (13)
% This part used enlarged CNT unit cell with 2*Nu carbon atoms in it.

% v_FT=zeros(nk,2*Nu,2*Nu);
% 
% 
% for u=-floor(2*nr):floor(2*nr)
%     u
%     for s1=1:Nu
%         for s2=1:Nu
%             deltaR=u*t_vec+posA(s1,:)-posA(s2,:);
%             v_FT(:,s1,s2)=v_FT(:,s1,s2)+(exp(1i*(q_vec*(K2*transpose(deltaR))+mu*(K1*transpose(deltaR))))* ...
%                 Upp/(sqrt((4*pi*eps0/q0^2*Upp*norm(deltaR))^2+1)))';
%             deltaR=u*t_vec+posA(s1,:)-posB(s2,:);
%             v_FT(:,s1,Nu+s2)=v_FT(:,s1,Nu+s2)+(exp(1i*(q_vec*(K2*transpose(deltaR))+mu*(K1*transpose(deltaR))))* ...
%                 Upp/(sqrt((4*pi*eps0/q0^2*Upp*norm(deltaR))^2+1)))';
%             deltaR=u*t_vec+posB(s1,:)-posA(s2,:);
%             v_FT(:,Nu+s1,s2)=v_FT(:,Nu+s1,s2)+(exp(1i*(q_vec*(K2*transpose(deltaR))+mu*(K1*transpose(deltaR))))* ...
%                 Upp/(sqrt((4*pi*eps0/q0^2*Upp*norm(deltaR))^2+1)))';
%             deltaR=u*t_vec+posB(s1,:)-posB(s2,:);
%             v_FT(:,Nu+s1,Nu+s2)=v_FT(:,Nu+s1,Nu+s2)+(exp(1i*(q_vec*(K2*transpose(deltaR))+mu*(K1*transpose(deltaR))))* ...
%                 Upp/(sqrt((4*pi*eps0/q0^2*Upp*norm(deltaR))^2+1)))';
%         end;
%     end;
% end;
% 
% v_FT=v_FT/nk;
% for i=-(nk-1)/2:(nk-1)/2
%     v_q(i+(nk+1)/2)=1/(4*Nu^2)*sum(sum(v_FT(i+(nk+1)/2,:,:)));
% end;
% 
% eps_q=1+v_q.*PI;
% fig=fig+1; figure(fig); hold on; box on;
% plot(q_vec,real(eps_q),'-k','LineWidth',3);
% ylabel('eps_q');
% xlim([0,4e9]);