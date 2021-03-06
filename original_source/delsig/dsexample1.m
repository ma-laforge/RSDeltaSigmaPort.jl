function dsm = dsexample1(dsm, LiveDemo)
%dsm = dsexample1(dsm, LiveDemo=0) Example lowpass/bandpass real/quadrature modulator design.
%Input dsm struct:
%   order
%   M
%   osr
%   f0
%   form
%   quadrature
% If quadrature == 0
%   Hinf
%   opt
% If quadrature == 1
%   NG
%   ING
%
% Output dsm: Above fields plus
%   nlev
%   ntf
%   ABCD
%   umax
%   amp, snr
%   peak_snr
%   coefficients.a
%   coefficients.g
%   coefficients.b
%   coefficients.c

%Handle arguments and defaults
if nargin<2
    LiveDemo = 0;
    if nargin<1
        dsm = [];
    end
end
if ~isfield(dsm, 'quadrature') || dsm.quadrature == 0
    FieldsAndDefaults = {
        'order' 3
        'osr'   16
        'M'     8
        'f0'    0
        'Hinf'  2
        'form'  'CRFB'
        'opt'   1
        'Atest' -3
        'Ftest' []
        'quadrature'   0
        };
else
    FieldsAndDefaults = {
        'order'	4
        'osr'	32
        'M'     8
        'f0'    1/16
        'form'  'PFB'
        'NG'    -50
        'ING'   -10
        'Atest' -3
        'Ftest' []
        'quadrature'   1
        };
end
for i = 1:size(FieldsAndDefaults,1)
    parameter = FieldsAndDefaults{i,1};
    if ~isfield(dsm, parameter)
        if isnan(FieldsAndDefaults{i,2})
            error('%s: The dsm argument must contain a field named ''%s''.', ...
                mfilename, parameter)
        else
            dsm.(parameter) = FieldsAndDefaults{i,2};
        end
    end
    eval([parameter '=dsm.(parameter);'])
end
% Computed defaults
if isempty(Ftest)
    if f0 == 0 || quadrature
       Ftest = 0.15/osr;
    else
       Ftest = f0 + 0.08/osr;
    end
end
% Derived parameters
nlev = M + 1;

clc
close all
if f0 == 0
    type = 'Lowpass';
else
    type = 'Bandpass';
end
if quadrature
    type = ['Quadrature ' type];
end


fprintf( 1, '\t%s-Order %s Example... ', ds_orderString(order,1), type );
fig_pos = { ...
    [  10 820  620 350]
    [ 650 820  620 350]
    [  10 440  620 350]
    [ 650 440  400 350]
    [  15 425  650 340]
    };
if LiveDemo
    sizes.lw = 2;       % LineWidth
    sizes.ms = 6;       % MarkerSize
    sizes.fs = 12;      % FontSize
    sizes.fw = 'bold';  % FontWeight
else
    sizes.lw = 1;       % LineWidth
    sizes.ms = 5;       % MarkerSize
    sizes.fs = 12;      % FontSize
    sizes.fw = 'normal';  % FontWeight
end

% NTF synthesis and realization
% fprintf(1,'Doing NTF synthesis and realization... ');
design_step = 1;
if ~quadrature
    ntf = synthesizeNTF(order,osr,opt,Hinf,f0);		% Optimized zero placement
    [a,g,b,c] = realizeNTF(ntf,form);
    z0 = exp(2i*pi*f0);
    b = [abs(b(1)+b(2)*(1-z0)) zeros(1,length(b)-1)];	% Use a single feed-in for the input
    ABCD = stuffABCD(a,g,b,c,form);
else
    ntf = synthesizeQNTF(order,osr,f0,NG,ING);
    ABCD = realizeQNTF(ntf,form,1);
end
% fprintf(1,'Done.\n');
figure(design_step); clf

set(design_step,'position',fig_pos{design_step});
DocumentNTF(ABCD,osr,f0,quadrature,sizes,1);
legend('NTF','STF', 'Location','SouthEast');
if LiveDemo
    drawnow
    pause
end

% Time-domain simulations
%fprintf(1,'Doing time-domain simulations... ');
design_step = design_step+1;
figure(design_step); clf; subplot('position', [0.07 0.15 0.88 0.75]);
set(gcf,'name','Time Domain')
set(gcf,'numbertitle','off');
set(gcf,'MenuBar','none');
set(design_step,'position',fig_pos{design_step});
% Time-domain plot
N = 100;
t = 0:N-1;
if ~quadrature
    u = undbv(Atest)*M*sin( 2*pi*Ftest*t );
    v = simulateDSM(u,ABCD,nlev);
    h = stairs( t, u, 'g');
    set(h,'Linewidth',sizes.lw);
    hold on;
    h = stairs( t, v, 'b');
    set(h,'Linewidth',sizes.lw);
    figureMagic([0 N],10,5, [-M M],2,1);
    axis([0 N -(M+0.25) M+0.25])
else
    u = undbv(Atest)*M*exp( 2i*pi*Ftest*t );
    v = simulateQDSM(u,ABCD,nlev);
    subplot(211);
    h = stairs( t, real(u), 'g');
    set(h,'Linewidth',sizes.lw);
    hold on;
    h = stairs( t, real(v), 'b');
    set(h,'Linewidth',sizes.lw);
    figureMagic([0 N],10,5, [-M M],2,1);
    subplot(212);
    h = stairs( t, imag(u), 'r');
    set(h,'Linewidth',sizes.lw);
    hold on;
    h = stairs( t, imag(v), 'm');
    set(h,'Linewidth',sizes.lw);
    figureMagic([0 N],10,5, [-M M],2,1);
end
xlabel('Sample Number');
if LiveDemo
    drawnow
    pause
end

% Example spectrum
design_step = design_step+1;
figure(design_step); clf; subplot('position', [0.13 0.15 0.82 0.75]);
set(gcf,'name','Example Spectrum')
set(gcf,'numbertitle','off');
set(gcf,'MenuBar','none');
set(design_step,'position',fig_pos{design_step});
PlotExampleSpectrum(ntf,M,osr,f0,quadrature,sizes,Atest,Ftest);
if LiveDemo
    drawnow
    pause
end

%figureMagic([0 0.5], 1/16, 2, [-140 0], 10, 2, [5.5 3], 'Example Spectrum');
% SQNR plot
design_step = design_step+1;
figure(design_step); clf;
set(design_step,'position',fig_pos{design_step});
if nlev==2
    [snr_pred,amp_pred] = predictSNR(ntf,osr);
    plot(amp_pred,snr_pred,'-');
    hold on;
end
[snr,amp] = simulateSNR(ABCD,osr,[],f0,nlev);
% fprintf(1,'Done.\n');
h = plot(amp,snr,'ob', 'Linewidth',sizes.lw);
set(h, 'MarkerSize',sizes.ms);
figureMagic([-120 0], 10, 2, [0 120], 10, 2, [3 3], 'SQNR vs. Input Level');
xlabel('Input Level (dBFS)');
ylabel('SQNR (dB)');
[peak_snr,peak_amp] = peakSNR(snr,amp);
msg = sprintf('Peak SQNR = %4.1fdB \n@ A = %4.1f dBFS',peak_snr,peak_amp);
h = text(-60,117,msg,'hor','cen', 'vertical','top');
set(h, 'FontSize',sizes.fs, 'FontWeight',sizes.fw);
% if LiveDemo
%     drawnow;
%     pause
% end

if quadrature    % Example I/Q mismatch
    fprintf(1,'Calculating effect of example I/Q mismatch... ');
    design_step = design_step+1;
    figure(design_step); clf;
    set(design_step,'position',fig_pos{design_step});
    ABCDr = mapQtoR(ABCD);
    ABCDr(2,end) = 1.01*ABCDr(2,end);   % 1% mismatch in first feedback
    [H G HI GI] = calculateQTF(ABCDr);
    f = ds_freq(osr,f0,quadrature);
    w = 2*pi*f;
    NTF  = squeeze( freqresp(H,w) );
    INTF = squeeze( freqresp(HI,w) );
    STF  = squeeze( freqresp(G,w) );
    ISTF = squeeze( freqresp(GI,w) );
    % fprintf(1,'Done.\n');
    plot(f,dbv(NTF),'b', 'Linewidth', 2);
    hold on;
    plot(f,dbv(INTF),'c', 'Linewidth', 2);
    plot(f,dbv(STF),'m', 'Linewidth', 2);
    plot(f,dbv(ISTF),'r', 'Linewidth', 2);
    inband = find( f>=f0-0.5/osr & f<=f0+0.5/osr );
    ng = dbp( mean(abs(NTF(inband)).^2) );
    plot(f0+0.5*[-1 1]/osr, ng*[1 1], 'k', 'Linewidth',3 );
    msg = sprintf('  NG = %.0fdB ', ng);
    h=text(f0+0.5/osr,ng,msg,'Hor','left','Vert','mid');
    set(h, 'FontSize',sizes.fs, 'FontWeight',sizes.fw);
    ing = dbp( mean(abs(INTF(inband)).^2) );
    plot(f0+0.5*[-1 1]/osr, ing*[1 1], 'k', 'Linewidth',3 );
    msg = sprintf('  rms(INTF) = %.0fdB ', ing);
    h=text(f0+0.5/osr,ing,msg,'Hor','left','Vert','mid');
    set(h, 'FontSize',sizes.fs, 'FontWeight',sizes.fw);
    irr = min( dbv(STF(inband)) - dbv(ISTF(inband)) );
    plot(f0+0.5*[-1 1]/osr, -irr*[1 1], 'k', 'Linewidth',3 );
    msg = sprintf('  IRR = %.0fdB ', irr);
    h=text(f0+0.5/osr,-irr,msg,'Hor','left','Vert','mid');
    set(h, 'FontSize',sizes.fs, 'FontWeight',sizes.fw);
    figureMagic([-0.5,0.5],1/16,2, [-80 15],10,2,[],'Example I/Q Mismatch');
    xlabel('frequency');
    legend('NTF','INTF','STF','ISTF',3);
    if LiveDemo
        drawnow
        pause
    end
end


if ~LiveDemo && ~quadrature   % Dynamic range scaling
    %fprintf(1,'Doing dynamic range scaling... ');
    design_step = design_step+1;
    figure(design_step); clf;
    set(design_step,'position',fig_pos{design_step});
    ABCD0 = ABCD;
    [ABCD, umax] = scaleABCD(ABCD0,nlev,f0);
    [a,g,b,c] = mapABCD(ABCD,form);
    % fprintf(1,'Done.\n');
    fprintf(1,'Verifying dynamic range scaling... ');
    u = linspace(0,0.95*umax,30);
    N = 1e4; N0 = 50;
    test_tone = cos(2*pi*f0*(0:N-1));
    test_tone(1:N0) = test_tone(1:N0) .* (0.5-0.5*cos(2*pi/N0*[0:N0-1]));
    maxima = zeros(order,length(u));
    for i = 1:length(u)
        ui = u(i);
        [v,xn,xmax] = simulateDSM( ui*test_tone, ABCD, nlev ); %#ok<ASGLU>
        maxima(:,i) = xmax(:);
        if any(xmax>1e2)
            fprintf(1,'Warning, umax from scaleABCD was too high.\n');
            umax = ui;
            u = u(1:i);
            maxima = maxima(:,1:i);
            break;
        end
    end
    % fprintf(1,'Done.\n');
    colors = hsv(order);
    cmd = 'legend( handles';
    handles = zeros(1,order);
    for i = 1:order
        handles(i) = plot(u,maxima(i,:),'--','color',colors(i,:));
        if i==1
            hold on;
        end
        cmd = [ cmd ',''State ' num2str(i) '''' ]; %#ok<AGROW>
        plot(u,maxima(i,:),'o','color',colors(i,:));
    end
    hold off; grid on;
    text(umax/2,0.05,'DC Input','Hor','cen','Vert','mid');
    figureMagic([ 0 umax],[],[], [0 1],0.1,2, [3 3],'State Maxima');
    set(gca,'position', [0.1 0.07, 0.85, 0.85]);
    cmd = [cmd ',4);'];
    eval(cmd);
else
    umax = [];
end


% The next step would be to size capacitors for each integrator state based
% on noise (kT/C) limits and area.

% Simulations modelling the effects of finite op-amp gain and capacitor
% errors should also be performed.

dsm.nlev = nlev;
dsm.ntf = ntf;
dsm.ABCD = ABCD;
dsm.amp = amp;
dsm.snr = snr;
dsm.umax = umax;
dsm.peak_snr = peak_snr;
if ~LiveDemo && ~quadrature
    dsm.coefficients.a = a;
    dsm.coefficients.g = g;
    dsm.coefficients.b = b;
    dsm.coefficients.c = c;
end
fprintf(1,'Done.\n');




