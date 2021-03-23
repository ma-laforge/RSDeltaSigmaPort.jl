%% Lowpass and Bandpass Demonstration Script
% Type cmd-shift-return to evaluate each cell in turn
close all; clear
LiveDemo = 1;

%% 2nd-order lowpass
dsm.order = 2;
dsm.osr = 16;
dsm.opt = 0;
dsm.Hinf = 2;
dsm.M = 8;      % nlev = M+1
dsexample1(dsm, LiveDemo); 

%% 5th-order lowpass
dsm.order = 5;
dsexample1(dsm, LiveDemo); 

%% 5th-order lowpass with optimized zeros
dsm.opt = 1;
dsexample1(dsm, LiveDemo); 

%% 5th-order lowpass with optimized zeros and larger Hinf
dsm.Hinf = 3;
dsexample1(dsm, LiveDemo); 

%% 7th-order lowpass; Hinf=2
dsm.order = 7;
dsm.osr = 8;
dsm.opt = 1;
dsm.Hinf = 2;
dsm.M = 16;      % nlev = M+1
dsm.Atest = -6;
dsexample1(dsm, LiveDemo);

%% 7th-order lowpass; Hinf=8
dsm.Hinf = 8;
dsexample1(dsm, LiveDemo); 

%% 6th-order bandpass
clear
LiveDemo = 1;
dsm.order = 6;
dsm.osr = 16;
dsm.opt = 1;
dsm.Hinf = 2;
dsm.M = 8;      
dsm.f0 = 1/6;   % Normalized center frequency
dsexample1(dsm, LiveDemo); 


