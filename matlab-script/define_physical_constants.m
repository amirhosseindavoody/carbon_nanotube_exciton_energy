%% Define physical constants
eV=1.6e-19; %[Joules]
hb=6.5e-16; %[eV.s]
aCC=1.42e-10; %[meters] carbon-carbon separation in graphene sheet.
aL=sqrt(3)*aCC; %lattice constant in graphene sheet.

e2p=0; %[Joules] parameter for band structure calculation
% t0=-3.033*eV; %[Joules] parameter for band structure calculation
t0=2.7*eV; %[Joules] parameter for band structure calculation
% s0=0.129; % parameter for band structure calculation
s0=0; % debug perposes

Zeff=4.14; %effective nuclear charge of carbon based on Nishant's paper
a_bohr=5.29e-11; %[meters] Bohr radius

Upp=11.3*eV; %the energy cost to place two electrons on a single site
eps0=8.85e-12; %permittivity of free space
q0=1.6e-19; %charge of electron

kappa=2; %dielectric constant due to core electrons in CNT

E_th=E_th*eV;