function [E_v,E_c,C_v,C_c]=graphene_energy(e2p,t0,s0,a1,a2,k)
f_k=exp(1i*k*transpose((a1+a2)/3))+exp(1i*k*transpose((a1-2*a2)/3))+exp(1i*k*transpose((a2-2*a1)/3));
% H=[e2p,t0*f_k;t0*conj(f_k),e2p];
% S=[1,s0*f_k;s0*conj(f_k),1];
% [C,E]=eig(H,S);
% E_v=E(1,1);
% E_c=E(2,2);
% C_v=C(:,1);
% C_c=C(:,2);
% C_v=C_v*(conj(C_v(1))/sqrt(conj(C_v(1))*C_v(1)));
% C_c=C_c*(conj(C_c(1))/sqrt(conj(C_c(1))*C_c(1)));

E_c=+t0*abs(f_k);
C_c(1,1)=+1/sqrt(2);
C_c(2,1)=+1/sqrt(2)*conj(f_k)/abs(f_k);

E_v=-t0*abs(f_k);
C_v(1,1)=+1/sqrt(2);
C_v(2,1)=-1/sqrt(2)*conj(f_k)/abs(f_k);

end

