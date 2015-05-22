%% This file reads dielectric function and Fourier transform of
% the Coulomb interaction based on the equations (13) and (17) in PRB 75,035407 (2007)

for iq=iq_min:iq_max
    q_vec(iq-iq_min+1)=iq*dk;
end;

load(filename,'v_FT','eps_q','v_q','ik_min','ik_max','mu');