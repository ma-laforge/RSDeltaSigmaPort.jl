%function adc = dsexample2
% Design example for a continuous-time lowpass DS ADC

% Design parameters
order = 3;
osr = 100;
nlev = 2;
f0 = 0;
opt = 0;
Hinf = 1.3;
umax = 0.83;
% Parameters for the continuous-time implementation
tdac = [0 1];	% DAC timing. [0 1] means zero-delay non-return-to-zero
form = 'FB';

% Derived parameters
M = nlev-1;

close all; clc
fprintf(1,'\t\t\t%dth-Order Continuous-Time Lowpass Example\n\n',order);
fig_pos = { [  10 595  600 215],
            [ 640 595  600 215],
            [  10 345  300 205],
            };

% NTF synthesis and realization
fprintf(1,'Doing NTF synthesis... ');
design_step = 1;
ntf0 = synthesizeNTF(order,osr,opt,Hinf,f0);		% Optimized zero placement
figure(design_step); clf
set(design_step,'position',fig_pos{design_step});
ntf_axes = DocumentNTF(ntf0,osr,f0);
drawnow;
%fprintf(1,'Paused... '); pause
fprintf(1,'Done.\n');

% Time-domain simulations
fprintf(1,'Doing time-domain simulations... ');
design_step = design_step+1;
figure(design_step); clf;
set(design_step,'position',fig_pos{design_step});
% Example spectrum
subplot('position', [0.05,0.1,0.6 0.8]);
Ftest = 0.3 * 0.5/osr /3;  % ~1/3 of the way across the passband
PlotExampleSpectrum(ntf0,M,osr,f0,[],[],[],Ftest,2^16);
title('Example Spectrum');
% SQNR plot
subplot('position',[.74 .18 .25 .65]);
if nlev==2
    [snr_pred,amp_pred] = predictSNR(ntf0,osr);
    plot(amp_pred,snr_pred,'-');
    hold on;
end
[snr,amp] = simulateSNR(ntf0,osr,[],f0,nlev);
plot(amp,snr,'og');
figureMagic([-100 0], 10, 2, [0 100], 10, 2, [7 3], 'Discrete-Time Simulation');
xlabel('Input Level (dBFS)');
ylabel('SQNR (dB)');
[peak_snr,peak_amp] = peakSNR(snr,amp);
msg = sprintf('peak SQNR = %4.1fdB  \n@ amp=%4.1fdB  ',peak_snr,peak_amp);
text(peak_amp-10,peak_snr,msg,'hor','right', 'vertical','middle');
msg = sprintf('OSR=%d ',osr);
text(0,5,msg,'hor','right');
title('SQNR Plot');
drawnow;
%fprintf(1,'Paused... '); pause
fprintf(1,'Done.\n');

% Continuous-Time Mapping
fprintf(1,'Mapping  to continuous-time... ');
design_step = design_step+1;
[ABCDc,tdac2] = realizeNTF_ct( ntf0, form, tdac);
[Ac Bc Cc Dc] = partitionABCD(ABCDc);
sys_c = ss(Ac,Bc,Cc,Dc);
%fprintf(1,'Paused... '); pause
fprintf(1,'Done.\n');
% Verify that the sampled pulse response of the CT loop filter
% matches the impulse response of the DT prototype
figure(design_step); clf;
set(design_step,'position',fig_pos{design_step});
set(design_step,'name','Continuous-Time Mapping')
set(design_step,'numbertitle','off');
set(design_step,'MenuBar','none');
n_imp = 10;
y = -impL1(ntf0,n_imp); % Negate impL1 to make the plot prettier
lollipop(0:n_imp,y);
hold on;
yl = floor(0.9*max(y)); plot([0.5 2], [yl yl], 'b', [0.5 2], [yl yl], 'bo');
text(2,yl,'   discrete-time');
grid on;
dt = 1/16;
yy = -pulse(sys_c,[0 0;tdac],dt,n_imp);
t = 0:dt:n_imp;
plot(t,yy,'g');
yl = floor(0.7*max(y)); plot([0.5 2], [yl yl], 'g');
text(2,yl,'  continuous-time');
title('Loop filter pulse/impulse responses (negated)');
% Map the cts system to its discrete time equivalent and check the NTF
sys_d = mapCtoD(sys_c,tdac);
ABCD = [sys_d.a sys_d.b; sys_d.c sys_d.d];
[ntf G] = calculateTF(ABCD);
ntf = cancelPZ(ntf);
axes(ntf_axes(1));
hold on;
plotPZ(ntf,'c',6);
hold off;
% Also plot the STF'
LF = zpk(sys_c);
L0 = LF(1);
f = linspace(0,0.5);
G = evalTFP(L0,ntf,f);
axes(ntf_axes(2));
hold on; 
plot(f, dbv(G), 'm');
hold off;
set(gcf,'name','NTF and STF')
if 0
    subplot('position',[.74 .25 .25 .5]);
    [f1x,f2x] = ds_f1f2(osr/1.5,f0);
    f = linspace(f1x,f2x);
    G = evalTFP(L0,ntf,f);
    plot(f*Fs/1e6, dbv(G), 'm');
    hold on; plot([f1 f2]*Fs/1e6, [-0.5 -0.5], 'k', 'LineWidth',3);
    figureMagic([f1x f2x]*Fs/1e6,0.5,2, [-0.5 0.5],.1,5, [6 2.5], 'NTF and STF Frequency Response');
    ylabel('dB')
    xlabel('MHz');
    title('Passband')
end

if 0
    design_step = design_step+1;
    fprintf(1,'Re-evaluating the SNR... ');
    figure(design_step); clf;
    set(design_step,'position',[830 15 400 400]);
    [snr amp] = simulateSNR(ABCD,osr,[],f0,nlev);
    fprintf(1,'Done.\n');
    plot(amp,snr,'o',amp,snr,'-')
    [peak_snr peak_amp] = peakSNR(snr,amp);
    msg = sprintf('Peak SNR \\approx %.0fdB at amp \\approx %-.0fdB',peak_snr,peak_amp);
    text(peak_amp-10, peak_snr, msg, 'hor', 'right');
    figureMagic([-100 0],10,1, [0 100],10,1, [3 3], 'SQNR vs. Input Amplitude');
    set(gca,'position',[0.12 0.15 0.8 0.7]);
    xlabel('Input Amplitude (dBFS)');
    ylabel('SNR (dB)');
    title('Continuous-Time Implementation');
end

design_step = design_step+1;
fprintf(1,'Doing dynamic range scaling... ');
% !!! This code assumes that the scale factors for the DT equivalent apply 
% !!! to the CT system. A system with an RZ DAC will have inter-sample peaks 
% !!! that exceed the values an the sampling instants.
[ABCDs,umax,S] = scaleABCD(ABCD,nlev,f0,1,[],umax,1e4);
S = S(1:order,1:order);	% Don't worry about the extra states used in the d-t model
Sinv = inv(S);
Acs=S*Ac*Sinv; Bcs=S*Bc;  Ccs=Cc*Sinv;
ABCDcs = [Acs Bcs; Ccs Dc];
sys_cs = ss(Acs,Bcs,Ccs,Dc);
%fprintf(1,'Paused... '); pause
fprintf(1,'Done.\n');
% ABCDcs needs to be checked with CT simulations to
% 1. Verify pulse response
% 2. Verify signal swings

adc.order = order;
adc.osr = osr;
adc.opt = opt;
adc.M = M;
adc.f0 = f0;
adc.ntf = ntf;
adc.ABCD = ABCD;
adc.umax = umax;
adc.peak_snr = peak_snr;
adc.form = form;
adc.ABCDc = ABCDc;
adc.L0  = L0;
adc.sys_c = sys_c;
adc.ABCDcs = ABCDcs;
adc.sys_cs = sys_cs;












