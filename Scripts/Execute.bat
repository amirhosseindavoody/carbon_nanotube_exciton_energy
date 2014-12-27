@echo off
title Run exciton code
cls

REM use the following variable to change back to the starting directory
set curr_dir=%cd%

set n_ch= 9
set m_ch= 7
set flg_dielectric= 1
set nkg= 501
set E_th= 1.5

REM change working directory to data storing directory
chdir /d C:\Users\Amirhossein\Google Drive\Research\Exciton\Data

echo Chirality is (%n_ch%,%m_ch%)
echo Dielectric flag is %flg_dielectric%
echo nkg = %nkg%
echo E_th= %E_th%
start C:\Amirhossein\CarbonNanotube-ExcitonEnergy\Release\CNT_Exciton_Energy.exe ch %n_ch% %m_ch% flg_dielectric %flg_dielectric% nkg %nkg% E_th %E_th%

chdir /d %curr_dir%
pause