#RSDeltaSigmaPort SNR calculations
#-------------------------------------------------------------------------------

"""`snr = calculateSNR(hwfft,f,nsig=1)`

Estimate the signal-to-noise ratio, given the in-band bins of a (Hann-windowed)
fft and the location of the input signal (f>0).

 - For nsig=1, the input tone is contained in `hwfft(f:f+2)`;
   this range is appropriate for a Hann-windowed fft.
 - Each increment in `nsig` adds a bin to either side.
 - The SNR is expressed in dB.
"""
function calculateSNR(hwfft, f::Int, nsig::Int = 1)
	signalBins = f-nsig+1:f+nsig+1
	signalBins = signalBins[signalBins .> 0]
	signalBins = signalBins[signalBins .<= length(hwfft)]
	s = norm(hwfft[signalBins]) #*4/(N*sqrt(3)) for true rms value
	noiseBins = collect(1:length(hwfft))
	deleteat!(noiseBins, signalBins)
	n = norm(hwfft[noiseBins])
	snr = (n!=0) ? dbv(s/n) : Inf
	return snr
end


"""`(peak_snr,peak_amp) = peakSNR(snr,amp)`

Estimate the snr peak by threading a line of slope 1 through the (amp,snr) data.
Both amp and snr are expressed in dB.
"""
function peakSNR(snr, amp)
	snr = deepcopy(snr); amp = deepcopy(amp)
	#Delete garbage data
	i = abs.(snr) .== Inf;   deleteat!(snr, i); deleteat!(amp, i)
	i = isnan.(snr);         deleteat!(snr, i); deleteat!(amp, i)
	i = snr .< 3;            deleteat!(snr, i); deleteat!(amp, i)
	n = length(amp)

	#Sort by amplitude
	i = sortperm(amp)
	amp = amp[i]; snr = snr[i]

	i = [true]
	local m = 0.0
	while any(i) && n > 3
		#Draw a 45-degree median line through the data
		tmp = sort(snr .- amp)
		if mod(n,2) == 0
			m = mean( tmp[div(n,2) .+ [1, 0]] )
		else
			m = tmp[div(n+1,2)]
		end
		#Discard data that is more than 6dB off 
		i = abs.(amp .- snr .+ m) .> 6
		deleteat!(snr, i); deleteat!(amp, i)
		n = length(amp)
	end

	peak_amp = maximum(amp)
	peak_snr = peak_amp + m
	return (peak_snr, peak_amp)
end


"""`(snr,amp,k0,k1,sigma_e2) = predictSNR(ntf,R=64,amp=...,f0=0)`

Predict the SNR curve of a binary delta-sigma modulator by using
the describing function method of Ardalan and Paulos.

The modulator is specified by a noise transfer function (ntf).
The band of interest is defined by the oversampling ratio (R)
and the center frequency (f0).
The input signal is characterized by the amp vector.
amp defaults to [-120 -110...-20 -15 -10 -9 -8 ... 0]dB, where 0 dB means
a full-scale (peak value = 1) sine wave. 

The algorithm assumes that the amp vector is sorted in increasing order;
once instability is detected, the remaining SNR values are set to -Inf.

# Output:
 - snr: a vector of SNR values (in dB)
 - amp: a vector of amplitudes (in dB)
 - k0: the quantizer signal gain
 - k1: the quantizer noise gain
 - sigma_e2: the power of the quantizer noise (not in dB)

The describing function method of A&P assumes that the quantizer processes
signal and noise components separately. The quantizer is modelled as two
(not necessarily equal) linear gains, k0 and k1, and an additive white
gaussian noise source of power sigma_e2. k0, k1 and sigma_e2 are calculated
as functions of the input.
The modulator's loop filter is assumed to have nearly infinite gain at
the test frequency.

Future versions may accommodate STFs.
"""
function predictSNR(ntf, R=64, amp=vcat(-120:10:-20, -15, -10:0); f0::Float64=0.0)
	Nb = 100
	XTAB = range(-2, stop=0, length=21)
	#The following code was used to create the YTAB matrix
	#YTAB = [];
	#for xi=XTAB
	#	YTAB = [YTAB; hypergeo(0.5,1,xi) hypergeo(0.5,2,xi)];
	#end
	YTAB = [0.46575960516930   0.67366999387741;
		0.47904652357101   0.68426650762558;
		0.49316295981407   0.69527947902679;
		0.50817364454269   0.70673173666000;
		0.52414894104004   0.71864765882492;
		0.54116523265839   0.73105299472809;
		0.55930554866791   0.74397552013397;
		0.57866013050079   0.75744456052780;
		0.59932720661163   0.77149158716202;
		0.62141352891922   0.78615015745163;
		0.64503526687622   0.80145609378815;
		0.67031890153885   0.81744754314423;
		0.69740217924118   0.83416539430618;
		0.72643494606018   0.85165339708328;
		0.75758063793182   0.86995816230774;
		0.79101717472076   0.88912981748581;
		0.82693856954575   0.90922164916992;
		0.86555624008179   0.93029111623764;
		0.90710091590881   0.95239937305450;
		0.95182400941849   0.97561222314835;
		1.00000000000000   1.00000000000000;
	]

	if f0==0
		band_of_interest = range(0, stop=pi/R, length=Nb)
	else
		band_of_interest = range(2*pi*(f0-0.25/R), stop=2*pi*(f0+0.25/R), length=Nb)
	end

	num, den = _zp2tf(ntf.z, ntf.p, 1)
	num1 = num - den

	N = length(amp)
	snr = zeros(1,N) .- Inf
	k0 = zeros(1,N)
	k1 = zeros(1,N)
	sigma_e2 = zeros(1,N)

	u = 10 .^ (amp/20)

	Nimp = 100
	unstable = false
	f_prev = 0.0; fprime = 0.0
	m0 = 0.0; rho = 0.0; erfinvu = 0.0
	for n=1:N
		#Calculate sigma_e2
		if f0==0
			erfinvu = erfinv(u[n])
			sigma_e2[n] = 1 - u[n]^2 - 2/pi*exp(-2*erfinvu^2)
		else #Sinusoidal input.
			#Solve sqrt(pi)*u/2 = rho * hypergeo(0.5,2,-rho^2);
			#Formulate as solve f(rho)=0, f = rho*M(0.5,2,-rho^2)-K
			#and use the secant method.
			K = 0.5 * sqrt(pi) * u[n]
			if n==1
				rho = u[n]^2 #Initial guess; otherwise use previous value.
				fprime = 1.0
			end
			drho = 1.0
			for itn = 1:20
				m0 = interp1_cubic(XTAB,YTAB[:,2],-rho^2)
				f = rho*m0 - K
				if( itn >1 )
					fprime = max((f-f_prev)/drho,0.5) #Secant approx.
				end
				if abs(f) < 1e-8; break; end #!Converged
				drho = -f/fprime
				if abs(drho) > 0.2; drho = sign(drho)*0.2; end
				if abs(drho) < 1e-6; break; end #!Converged
				rho = rho + drho
				f_prev = f
			end
			m1 = interp1_cubic(XTAB,YTAB[:,1], -rho^2)
			sigma_e2[n] = 1 - u[n]^2/2 - 2/pi*m1^2
		end #Sinusoidal input

		#Iterate to solve for k1 and sigma_1.
		#Using one of MATLAB's nonlinear equation solvers would be more efficient,
		#but this function code would then require the optimization toolbox.
		#!Future work: put in 2-D BFGS code.
		k1[n] = (n>1) ? k1[n-1] : 1.2 #Use the previous value of k1 as the initial guess.
		k1_prev = 0
		itn = 0
		k1sigma1 = (f0==0) ? sqrt(2/pi)*exp(-erfinvu^2) : sqrt(2/pi)*m1
		local sigma_1 = 0.0
		while abs(k1[n]-k1_prev) > 1e-6*(1+k1[n]) && itn < 100
			#Create the function: H_hat = L1/(1-k1*L1)=(H-1)/(H*(1-k1)+k1).
			den1 = (1-k1[n])*num + den*k1[n]
			#Calculate pGain, the square of the 2-norm of H_hat.
			pGain, Nimp = powerGain(num1, den1, Nimp)
			if isinf(pGain)
				unstable = true
				break
			end

			sigma_1 = sqrt(pGain*sigma_e2[n])
			k1_prev = k1[n]
			k1[n] = k1sigma1/sigma_1
			itn = itn+1
		end
		if unstable
			break
		end

		if f0==0 
			y0 = sqrt(2)*erfinvu*sigma_1
			k0[n] = u[n]/y0
		else
			k0[n] = sqrt(2/pi)*m0/sigma_1
		end

		h = _freqz(num, (1-k1[n])*num + k1[n]*den, band_of_interest)
		#For both DC and sine wave inputs, use u^2/2 as the signal 
		#power since true DC measurements are usually impossible.
		snr[n] = dbp( 0.5*u[n]^2 / (sum(h.^2)/(R*Nb)*sigma_e2[n]) )
	end

	return (snr, amp, k0, k1, sigma_e2)
end

"""`(snr,amp) = simulateSNR(ΔΣmodel,osr,amp; f0=0,nlev=2,f=1/(4*osr),k=13,quadrature=false)`

Determine the SNR for a delta-sigma modulator by using simulations.

`ΔΣmodel` describes the modulator by a noise transfer function (ntf) and the
number of quantizer levels (`nlev`). Alternatively, `ΔΣmodel` may be a
function taking the input signal as its sole argument:
 - `ΔΣmodel::ZPKData`: z, p, k transfer function
 - `ΔΣmodel::Array` ABCD matrix
 - `ΔΣmodel::Function` funΔΣ(u) generating ΔΣ output.

The band of interest is defined by the oversampling ratio (`osr`)
and the center frequency (`f0`).
The input signal is characterized by the `amp` vector and the `f` variable.

 - `amp` defaults to [-120 -110...-20 -15 -10 -9 -8 ... 0]dB, where 0 dB means
   a full-scale (peak value = nlev-1) sine wave.
 - `f` is the input frequency, normalized such that 1 -> fs;
 - `f` is rounded to an FFT bin.

Using sine waves located in FFT bins, the SNR is calculated as the ratio
of the sine wave power to the power in all in-band bins other than those
associated with the input tone. Due to spectral smearing, the input tone
is not allowed to lie in bins 0 or 1. The length of the FFT is 2^k.

 - If ntf is complex, simulateQDSM (which is slow) is called.
 - If ABCD is complex, simulateDSM is used with the real equivalent of ABCD
   in order to speed up simulations.

Future versions may accommodate STFs.
"""
function simulateSNR(arg1, osr::Int=64, amp=vcat(-120:10:-20, -15, -10:0);
	f0::Float64=0.0, nlev::Int=2, f::Float64=NaN, k::Int=13, quadrature::Bool=false
)
	is_function = isa(arg1, Function)

	#Look at arg1 and decide if the system is quadrature
	quadrature_ntf = false
	if !is_function
		if isa(arg1, ZPKData)
			pz = [arg1.p, arg1.z]
			for i=1:2
				if any( abs.(imag.(poly(pz[i]))) .> 1e-4 )
					quadrature = true
					quadrature_ntf = true
				end
			end
		else #ABCD matrix
			if !all(all(imag.(arg1) .== 0))
				quadrature = true
			end
		end
	end

	osr_mult = 2
	if f0!=0 && !quadrature
   	osr_mult = 2*osr_mult
	end
	if isnan(f)
		f = f0 + 0.5/(osr*osr_mult) #Halfway across the band
	end
	M = nlev-1
	if quadrature && !quadrature_ntf	
		#Modify arg1 (ABCD) and nlev so that simulateDSM can be used
		nlev = [nlev; nlev]
		arg1 = mapQtoR(arg1)
	end

	if abs(f-f0) > 1/(osr*osr_mult)
		@warn("The input tone is out-of-band.")
	end

	N = 2^k
	if N < 8*2*osr #Require at least 8 bins to be "in-band"
		@warn("Increasing k to accommodate a large oversampling ratio.")
		k = ceil(Int, log2(8*2*osr))
		N = 2^k
	end
	F = round(Int, f*N)
	if abs(F)<=1
		@warn("Increasing k to accommodate a low input frequency.")
		#Want f*N >= 1
		k = ceil(Int, log2(1/f))
		N = 2^k
		F = 2
	end

	Ntransient = 100
	soft_start = 0.5*(1 .- cos.(2π/Ntransient * (0:div(Ntransient,2)-1)))
	if !quadrature
		tone = M * sin.(2π*F/N * (0:(N+Ntransient-1)))
		tone[1:div(Ntransient,2)] = tone[1:div(Ntransient,2)] .* soft_start
	else
		tone = M * exp.(2π*j*F/N * (0:(N+Ntransient-1)) )
		tone[1:div(Ntransient,2)] = tone[1:div(Ntransient,2)] .* soft_start
		if !quadrature_ntf
			tone = [real(tone); imag(tone)]
		end
	end
	window = .5*(1 .- cos.(2π*(0:N-1)/N) ) #Hann window
	if f0==0
		#Exclude DC and its adjacent bin
		inBandBins = div(N,2) .+ (3:round(Int, N/(osr_mult*osr)))
		F = F-2
	else
		f1 = round(Int, N*(f0-1/(osr_mult*osr)))
		f2 = round(Int, N*(f0+1/(osr_mult*osr)))
		inBandBins = div(N,2) .+ (f1:f2) #Should exclude DC
		F = F-f1+1
	end

	snr = zeros(size(amp))
	i=1
	for A = 10 .^ (amp/20)
		if is_function
			(v,) = arg1(A*tone)
		elseif quadrature_ntf
			(v,) = simulateQDSM(A*tone, arg1, nlev)
		else
			(v,) = simulateDSM(A*tone, arg1, nlev)
			if quadrature
				v = v[1,:] + im*v[2,:]
			end
		end
		hwfft = fftshift(fft(window .* v[1+Ntransient:N+Ntransient]))
		snr[i] = calculateSNR(hwfft[inBandBins], F)
		i=i+1
	end

	return (snr,amp)
end

#Last line
