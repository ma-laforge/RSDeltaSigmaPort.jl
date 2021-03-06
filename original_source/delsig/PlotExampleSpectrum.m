function PlotExampleSpectrum(ntf,M,osr,f0,quadrature,sizes,Atest,Ftest,N)
%PlotExampleSpectrum(ntf|mod_struct,M=1,osr=64,f0=0,quadrature=0,sizes,Atest,Ftest,N)
% ntf|mod_struct is either the NTF or a struct containing 
% ntf, M=1, osr=64, f0=0, quadrature=0

% Handle the input arguments
parameters = {'ntf' 'M' 'osr' 'f0' 'quadrature' 'sizes' 'Atest' 'Ftest' 'N'};
defaults = { [] 1 64 0 0 NaN -3 [] 2^12};
for arg_ii=1:length(defaults)
    parameter = parameters{arg_ii};
    if arg_ii>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{arg_ii};'])
    end
end
if nargin==1 & isstruct(ntf)
    dsm = ntf;
    ntf = dsm.ntf;
    fields = {'M','f0','quadrature'};
    for i = 1:length(fields)
        field = fields{i};
        if isfield(dsm,field)
            eval([field '= dsm.' field ';'])
        end
    end
end
if isnumeric(sizes) & isnan(sizes)
    clear sizes;
    sizes.lw = 1;        % LineWidth
    sizes.ms = 5;        % MarkerSize
    sizes.fs = 12;       % FontSize
    sizes.fw = 'normal'; % FontWeight
end
[f1 f2] = ds_f1f2(osr,f0,quadrature);
% Computed defaults
if isempty(Ftest)
    if f0 == 0 || quadrature
       Ftest = 0.15/osr;
    else
       Ftest = f0 + 0.08/osr;
    end
end

delta = 2;

% Plot an example spectrum
Amp = undbv(Atest);  % Test tone amplitude, relative to full-scale.
f1_bin = round(f1*N);
f2_bin = round(f2*N);
fin = round(Ftest*N);
if ~quadrature
    u = Amp*M*cos((2*pi/N)*fin*[0:N-1]);
    v = simulateDSM(u,ntf,M+1);
else
    u = Amp*M*exp((2i*pi/N)*fin*[0:N-1]);
    v = simulateQDSM(u,ntf,M+1);
end
window = ds_hann(N);
NBW = 1.5/N;
spec0 = fft(v.*window)/(M*N/4);
if ~quadrature
    freq = linspace(0,0.5,N/2+1);
    if f0==0
        semilogx(freq,dbv(spec0(1:N/2+1)),'c','Linewidth',1);
    else
        plot(freq,dbv(spec0(1:N/2+1)),'c','Linewidth',1);
    end
	hold on
    % circular smoothing; not re-scaled
	spec_smoothed = circ_smooth(abs(spec0).^2, 16);
	plot(freq, dbp(spec_smoothed(1:N/2+1)), 'b', 'Linewidth', 3);
	Snn = abs(evalTF(ntf,exp(2i*pi*freq))).^2 * 2*2/12 * (delta/M)^2;
	plot(freq, dbp(Snn*NBW), 'm', 'Linewidth', 1);
    snr = calculateSNR(spec0(f1_bin+1:f2_bin+1),fin-f1_bin);
    msg = sprintf(' SQNR = %.1fdB\n @ A = %.1fdBFS\n  & OSR = %.0f\n', ...
      snr, dbv(spec0(fin+1)), osr );
    if f0<0.25
        h = text(f0+0.5/osr, -5, msg, 'Hor','Left', 'Ver','top');
    else
        h = text(f0-0.5/osr, -5, msg, 'Hor','Right', 'Ver','top');
    end
    set(h, 'FontSize',sizes.fs, 'FontWeight',sizes.fw);
    h = text(0.5,-135,sprintf('NBW=%.1e ',NBW),'hor','right');
    set(h, 'FontSize',0.8*sizes.fs, 'FontWeight',sizes.fw);
    figureMagic([0 0.5],1/16,4, [-140 0],10,2);
else
    spec0 = fftshift(spec0/2);
    freq = linspace(-0.5,0.5,N+1); freq(end)=[];
	plot(freq,dbv(spec0),'c','Linewidth',1);
	hold on
	spec_smoothed = circ_smooth(abs(spec0).^2, 16);
	plot(freq, dbp(spec_smoothed), 'b', 'Linewidth', 3);
	Snn = abs(evalTF(ntf,exp(2i*pi*freq))).^2 * 2/12 * (delta/M)^2;
	plot(freq, dbp(Snn*NBW), 'm', 'Linewidth', 1);
    snr = calculateSNR(spec0(N/2+1+[f1_bin:f2_bin]),fin-f1_bin);
    msg = sprintf('SQNR = %.1fdB\n @ A=%.1fdBFS & osr=%.0f\n', ...
      snr, dbv(spec0(N/2+fin+1)), osr );
    if f0>=0
        text(f0-0.05, -15, msg,'Hor','Right');
    else
        text(f0+0.05, -15, msg,'Hor','Left');
    end
    text(-0.5,-135,sprintf(' NBW=%.1e',NBW),'hor','left');
    figureMagic([-0.5 0.5],0.125,2, [-140 0],10,2);
end
xlabel('Normalized Frequency');
ylabel('PSD (dBFS/NBW)');
% printmif(fullfile('MIF','spec'), [5 2], 'Helvetica10')



