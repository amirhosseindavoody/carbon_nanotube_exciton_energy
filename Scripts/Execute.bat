@echo off
title Run exciton code
cls
set n_ch= 16
set m_ch= 0
set flg_dielectric= 0
set nkg= 501
set E_th= 1.5
echo Chirality is (%n_ch%,%m_ch%)
echo Dielectric flag is %flg_dielectric%
echo nkg = %nkg%
echo E_th= %E_th%
start CNT_Exciton_Energy.exe ch %n_ch% %m_ch% flg_dielectric %flg_dielectric% nkg %nkg% E_th %E_th%
pause