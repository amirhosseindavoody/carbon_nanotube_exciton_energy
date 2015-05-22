%% This code calculates the band structure of single wall carbon nanotubes through simple tight-binding method.
% Amirhossein Davoody
% Last modified: 1/9/2014

clear all;  clc; fig=0;
close all;
dir='C:\Users\Amirhossein\Documents\My Dropbox\Research\Exciton\data\';

%% Input simulation parameters
nC=6; %chiral vector coefficient
mC=5; %chiral vector coefficient 0<mC<nC

nkg=501; %number of mesh points in the first brillioun zone of graphene. [This number should be set to an odd number.]
nr=100; %number of mesh points in real space in each dimension
E_th=2; %threshold energy in eV units.
n_sub=2; %number of subbands considered in calculation of excitonic states
Kcm_max=1.5e9; %maximum amount of Kcm in units of nanometers.
flg_dielectric=0; %when set to 1 dielectric function is calculated and stored in file.


%% Initialize simulation
define_physical_constants
calculate_geometrical_properties

filename=[dir,'dielectric','_chir_(',num2str(nC),',',num2str(mC),')_nkg_',num2str(nkg),'_nr_',num2str(nr),...
    '_Eth_',num2str(E_th/eV),'_nsub_',num2str(n_sub),'_s0_',num2str(s0),'_t0_',num2str(t0/eV),'_Kcm-max_',num2str(Kcm_max*1e-9),'.mat'];

%% Calculate carbon nanotube band structure
calculate_cnt_bandstructure
return
%% Define boundaries of k-space
iKcm_min=-floor(Kcm_max/dk);
iKcm_max=+floor(Kcm_max/dk);
ikr_upper=(iKcm_max-ik_min)+iKcm_max;
ikr_lower=-ikr_upper;
tmp=iKcm_max-ik_min;
iq_min=-(2*tmp);
iq_max=+(2*tmp);

% return;

%% Calculate eps(mu,q)
if (flg_dielectric==1)
    calculate_dielectric_function
    mu=abs(min_sub(1));
    save(filename,'v_FT','eps_q','v_q','ik_min','ik_max','mu');
else
    load_dielectric_function
end;
nq=numel(q_vec);

% mu_vec= -Nu+1:Nu-1;
% fig=fig+1; figure(fig); hold on; box on;
% surf(q_vec,mu_vec,real(eps_q),'EdgeColor','none');
% axis tight;
% title('\epsilon(q)');

% fig=fig+1; figure(fig); hold on; box on;
% surf(q_vec,mu_vec,real(v_q)/eV,'EdgeColor','none');
% axis tight;
% title('v(q)');

% PI=(eps_q-1)./v_q;
% fig=fig+1; figure(fig); hold on; box on;
% surf(q_vec,mu_vec,real(PI)*eV,'EdgeColor','none');
% axis tight;
% title('\Pi(q)');

%%
% fig=fig+1; figure(fig); hold on; box on;
% for mu=0:0
%     mu
%     tmp=zeros(1,nq);
%     tmp(1,:)=(eps_q(mu+Nu,:));
%     plot(q_vec,abs(tmp),'-b','LineWidth',3);
% end;
% xlabel('q [1/m]')
% axis tight;

% return;
%% Calculate self-energy corrections
calculate_self_energies

% return;
%% Calculate wavefunction of Psi with Kcm=0
% mu_cm=0;
% iKcm=0;
% Kcm=mu_cm*K1+iKcm*dk*K2;
% calculate_kernel_half_matrix;
% Psi0_A2=Psi0_A2_tmp;
% plot_exciton_wavefunction
% 
% return;
%% Calculate Interaction Kernels
mu_cm=0;
Kcm_vec=dk*(iKcm_min:1:iKcm_max);
nK_cm=numel(Kcm_vec);
nkr_max=2*ikr_upper+1;

Ex_A1=zeros(nkr_max,nK_cm);
Ex0_A2=zeros(nkr_max,nK_cm);
Ex1_A2=zeros(nkr_max,nK_cm);
Ef=zeros(nkr_max,nK_cm);
Psi_A1=zeros(nkr_max,nkr_max,nK_cm);
Psi0_A2=zeros(nkr_max,nkr_max,nK_cm);
Psi1_A2=zeros(nkr_max,nkr_max,nK_cm);

tmp10=zeros(n_sub,n_sub,nK_cm,2,2);
tic
% for iKcm=iKcm_min:iKcm_max
for iKcm=iKcm_min:iKcm_max
    iKcm
    Kcm=mu_cm*K1+iKcm*dk*K2;
    calculate_kernel_half_matrix
    Ef(1:nk,iKcm-iKcm_min+1)=Ef_tmp(1:nk);
    Ex_A1(1:nkr,iKcm-iKcm_min+1)=Ex_A1_tmp(1:nkr);
    Ex0_A2(1:nkr,iKcm-iKcm_min+1)=Ex0_A2_tmp(1:nkr);
    Ex1_A2(1:nkr,iKcm-iKcm_min+1)=Ex1_A2_tmp(1:nkr);
    Psi_A1(1:nkr,1:nkr,iKcm-iKcm_min+1)=Psi_A1_tmp;
    Psi0_A2(1:nkr,1:nkr,iKcm-iKcm_min+1)=Psi0_A2_tmp;
    Psi1_A2(1:nkr,1:nkr,iKcm-iKcm_min+1)=Psi1_A2_tmp;
end;
filename=[dir,'Dispersion','_chir_(',num2str(nC),',',num2str(mC),')_nkg_',num2str(nkg),'_nr_',num2str(nr),...
    '_Eth_',num2str(E_th/eV),'_nsub_',num2str(n_sub),'_s0_',num2str(s0),'_t0_',num2str(t0/eV),'_Kcm-max_',num2str(Kcm_max*1e-9),'.mat'];
save(filename,'Ex_A1','Ex0_A2','Ex1_A2','Psi_A1','Psi0_A2','Psi1_A2');

toc

%% Plot the energy dispersion
tmp=min(Ex0_A2(1,:));
% tmp=0;

position3=[0.1,0.1,0.3,0.3];
position1=[0.1,0.5,0.3,0.3];
position2=[0.5,0.5,0.3,0.3];
position4=[0.5,0.1,0.3,0.3];

fig=fig+1; figure(fig);
% figure(4);
subplot(2,2,1); hold on; box on;
for i=1:nk
    if (min(Ex0_A2(i,:))<=min(Ef(1,:)))
        plot(2*Kcm_vec,(Ex0_A2(i,:)-tmp)/eV,'-b','LineWidth',3);
    end;
end;
axis tight;
xlabel('2K_{cm} [1/m]');
ylabel('E_{11}(A_2;S=0) [eV]');

subplot(2,2,2); hold on; box on;
for i=1:nk
    if (min(Ex_A1(i,:))<=min(Ef(1,:)))
        plot(2*Kcm_vec,(Ex_A1(i,:)-tmp)/eV,'-k','LineWidth',3);
    end;
end;
axis tight;
xlabel('2K_{cm} [1/m]');
ylabel('E_{11}(A_1;S=0,1) [eV]');

subplot(2,2,3); hold on; box on;
for i=1:nk
    if (min(Ex1_A2(i,:))<=min(Ef(1,:)))
        plot(2*Kcm_vec,(Ex1_A2(i,:)-tmp)/eV,'-r','LineWidth',3);
    end;
end;
axis tight;
xlabel('2K_{cm} [1/m]');
ylabel('E_{11}(A_2;S=1) [eV]');

%plot bottom of subbands
subplot(2,2,4); hold on; box on;
for i=1:nk
    if (min(Ex0_A2(i,:))<=min(Ef(1,:)))
        plot([0,1],[min(Ex0_A2(i,:)-tmp)/eV,min(Ex0_A2(i,:)-tmp)/eV],'-b','LineWidth',3);
    end;
end;

for i=1:nk
    if (min(Ex_A1(i,:))<=min(Ef(1,:)))
        plot([1.2,2.2],[min(Ex_A1(i,:)-tmp)/eV,min(Ex_A1(i,:)-tmp)/eV],'-k','LineWidth',3);
    end;
end;

for i=1:nk
    if (min(Ex1_A2(i,:))<=min(Ef(1,:)))
        plot([2.4,3.4],[min(Ex1_A2(i,:)-tmp)/eV,min(Ex1_A2(i,:)-tmp)/eV],'-r','LineWidth',3);
    end;
end;

% for i=1:nk
%     plot([3.6,4.6],[min(Ef(i,:)-tmp)/eV,min(Ef(i,:)-tmp)/eV],'-g','LineWidth',3);
% end;
axis tight;