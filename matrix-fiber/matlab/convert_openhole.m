

clear all
clc
% inputPath='..\CAE\notched.inp';
% inputPath='..\CAE\vuel_test.inp';
% inputPath='..\CAE\onecrack.inp';
% inputPath='..\CAE\doublephase-test.inp';
% inputPath='..\CAE\doublephase.inp';
% MatProp=[1190 26.5e9 2.6e9 1.3e9 0 0.35 0.43 20.25e6 17.28e6 17.28e6 3000e6 0.622E3 ...
%     0.472E3 0.472E3 62E3 3e-4 3e-4 0.06 0.06 1 1.732 0];
% MatProp=[1190 26.5e9 26.5e9 11e9 11e9 0.3 0.3 20.25e9 17.28e9 17.28e9 1000e9 0.622E9 ...
    % 0.472E9 0.472E9 30E3 3e-4 3e-4 0.03 0.03 1 1 0];
% Abaqus2PhasefieldUEL3D(inputPath,MatProp,"double")



%%%%%%%%%%%%%%%%%%%%%%%% open hole %%%%%%%%%%%%%%%%%%%%%%
inputPath='..\CAE\openhole.inp';
MatProp=[1190 161e9 11.4e9 5.17e9 3.98e9 0.32 0.43 60e6 50e6 90e6 2806e6 0.293E3 0.631E3 0.631E3 112.7E3 4e-4 4e-4 0.05 0.05 1 1 0;
    1190 161e9 11.4e9 5.17e9 3.98e9 0.32 0.43 60e6 50e6 90e6 2806e6 0.293E3 0.631E3 0.631E3 112.7E3 4e-4 4e-4 0.05 0.05 0 1 0;
    1190 161e9 11.4e9 5.17e9 3.98e9 0.32 0.43 60e6 50e6 90e6 2806e6 0.293E3 0.631E3 0.631E3 112.7E3 4e-4 4e-4 0.05 0.05 -1 1 0;
    1190 161e9 11.4e9 5.17e9 3.98e9 0.32 0.43 60e6 50e6 90e6 2806e6 0.293E3 0.631E3 0.631E3 112.7E3 4e-4 4e-4 0.05 0.05 -1 0 0;
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%main program
Abaqus2PhasefieldUEL3D(inputPath,MatProp,"double",true)