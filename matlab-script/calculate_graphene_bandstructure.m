%% Calculate graphene band structure
kx_vec=linspace(0,4*b1(1),nk);
ky_vec=linspace(0,4*b1(2),nk);

for i=1:nk
    for j=1:nk
        k=[kx_vec(i),ky_vec(j)];
        [E_v(i,j),E_c(i,j)]=graphene_energy(e2p,t0,s0,a1,a2,k);
    end;
end;
E_c=E_c/eV;
E_v=E_v/eV;
fig=fig+1; figure(fig); hold on;
surf(ky_vec,kx_vec,E_c,'EdgeColor','none');
axis equal; axis tight; box on;
set(gcf,'renderer','Zbuffer');
view(90,90);

fig=fig+1; figure(fig); hold on;
surf(ky_vec,kx_vec,E_v,'EdgeColor','none');
axis equal; axis tight; box on;
set(gcf,'renderer','Zbuffer');
view(90,90);