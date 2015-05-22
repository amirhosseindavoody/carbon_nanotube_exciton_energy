%% Define discretized reciprocal space of the carbon nanotube
dk=norm(b1)/(nkg-1);
ikc_max=+floor(pi/norm(t_vec)/dk);
ikc_min=-floor(pi/norm(t_vec)/dk);
k_vec=dk*(ikc_min:1:ikc_max); %k_vec is defined here for plotting reasons and will be over written at the end of this file.
nkc=numel(k_vec); %number of mesh points in the first Brillioun zone of the CNT.
K1=(-t2*b1+t1*b2)/Nu;
K2=(mC*b1-nC*b2)/Nu;
K2=K2/norm(K2);


%% Calculate carbon nanotube bandstructure
E_c=zeros(Nu,nkc);
E_v=zeros(Nu,nkc);
C_v=zeros(Nu,nkc,2);
C_c=zeros(Nu,nkc,2);
for mu=1-Nu/2:Nu/2
    for ik=-(nkc-1)/2:(nkc-1)/2
        k=mu*K1+ik*dk*K2;
        [E_v(mu+Nu/2,ik+(nkc+1)/2),E_c(mu+Nu/2,ik+(nkc+1)/2),C_v(mu+Nu/2,ik+(nkc+1)/2,:),C_c(mu+Nu/2,ik+(nkc+1)/2,:)]=graphene_energy(e2p,t0,s0,a1,a2,k);
    end;
end;

fig=fig+1; figure(fig); hold on;
plot(k_vec,E_v/eV,'-','LineWidth',3);
plot(k_vec,E_c/eV,'-','LineWidth',3);

%% This part finds the bands with minimum energy and the limits of k-vector the energy of that band is below threshold energy
tmp=min(E_c,[],2);
[i,j]=min(tmp);
min_sub(1)=-abs(j-Nu/2);
min_sub(2)=-min_sub(1);
[Ec_min(1),ik1]=min(E_c(min_sub(1)+Nu/2,:));
[Ec_min(2),ik2]=min(E_c(min_sub(2)+Nu/2,:));
ik1=ik1-(nkc+1)/2;
ik2=ik2-(nkc+1)/2;
clear tmp;

n=0;
m=0;
counter1=1;
flg_u=0;
flg_d=0;
while((flg_u==0)||(flg_d==0))
    if (flg_u==0)
        n=n+1;
        k=min_sub(1)*K1+(ik1+n)*dk*K2;
        [E_v_tmp,E_c_tmp,C_v_tmp,C_c_tmp]=graphene_energy(e2p,t0,s0,a1,a2,k);
        plot((ik1+n)*dk,E_c_tmp/eV,'r*');
        if ((E_c_tmp-Ec_min(1))>E_th)
            flg_u=1;
        else
            counter1=counter1+1;
            ik_max(1)=ik1+n;
        end;
    end;
    if (flg_d==0)
        m=m-1;
        k=min_sub(1)*K1+(ik1+m)*dk*K2;
        [E_v_tmp,E_c_tmp,C_v_tmp,C_c_tmp]=graphene_energy(e2p,t0,s0,a1,a2,k);
        plot((ik1+m)*dk,E_c_tmp/eV,'*');
        if ((E_c_tmp-Ec_min(1))>E_th)
            flg_d=1;
        else
            counter1=counter1+1;
            ik_min(1)=ik1+m;
        end;
    end;
end;

counter2=1;
n=0;
m=0;
flg_u=0;
flg_d=0;
while((flg_u==0)||(flg_d==0))
    if (flg_u==0)
        n=n+1;
        k=min_sub(2)*K1+(ik2+n)*dk*K2;
        [E_v_tmp,E_c_tmp,C_v_tmp,C_c_tmp]=graphene_energy(e2p,t0,s0,a1,a2,k);
        if ((E_c_tmp-Ec_min(2))>E_th)
            flg_u=1;
        else
            counter2=counter2+1;
            ik_max(2)=ik2+n;
        end;
    end;
    if (flg_d==0)
        m=m-1;
        k=min_sub(2)*K1+(ik2+m)*dk*K2;
        [E_v_tmp,E_c_tmp,C_v_tmp,C_c_tmp]=graphene_energy(e2p,t0,s0,a1,a2,k);
        if ((E_c_tmp-Ec_min(2))>E_th)
            flg_d=1;
        else
            counter2=counter2+1;
            ik_min(2)=ik2+m;
        end;
    end;
end;

if ((ik_max(1)-ik_min(1))~=(ik_max(2)-ik_min(2)))
        disp('Error in calculating lowest subbands!');
        exit;
end;

ik_min=min(ik_min(1),ik_min(2));
ik_max=-ik_min;
nk=ik_max-ik_min+1; %number of mesh points from the each subband considered in calculation of exciton energy.

k_vec=dk*(ik_min:1:ik_max);


%% Plot the bands that are calculated in the previous section
E_c=zeros(2,nk);
E_v=zeros(2,nk);
C_v=zeros(2,nk,2);
C_c=zeros(2,nk,2);

for ik=ik_min:ik_max
    mu=min_sub(1);
    k=mu*K1+ik*dk*K2;
    [E_v(1,ik-ik_min+1),E_c(1,ik-ik_min+1),C_v(1,ik-ik_min+1,:),C_c(1,-ik_min+1,:)]=graphene_energy(e2p,t0,s0,a1,a2,k);
    
    mu=min_sub(2);
    k=mu*K1+ik*dk*K2;
    [E_v(2,ik-ik_min+1),E_c(2,ik-ik_min+1),C_v(2,ik-ik_min+1,:),C_c(2,-ik_min+1,:)]=graphene_energy(e2p,t0,s0,a1,a2,k);
end;

fig=fig+1; figure(fig); hold on;
plot(k_vec,E_v(1,:)/eV,'-r','LineWidth',3);
plot(k_vec,flipud(E_c(1,:))/eV,'-b','LineWidth',3);
plot(k_vec,E_v(2,:)/eV,'--r','LineWidth',3);
plot(k_vec,flipud(E_c(2,:))/eV,'--b','LineWidth',3);
axis tight; box on;