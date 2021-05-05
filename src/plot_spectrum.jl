#RSDeltaSigmaPort: Generate frequency spectrum plots
#-------------------------------------------------------------------------------


#==Frequency-spectrum plots
===============================================================================#
function plotSpec()
	plot = cons(:plot, linlin, title = "Signal Spectrum", legend=true,
		xyaxes=set(xmin=0, xmax=0.5, ymin=-120, ymax=0),
		labels=set(xaxis="Normalized Frequency", yaxis="dBFS"),
	)
	return plot
end

function plotSpec!(plot, sig; M::Int=1, id::String="signal", color=:blue)
	sig = conv2seriesmatrix2D(sig)
	M, N = size(sig)
	𝑓pts = div(N,2)+1 #Only need half the spectrum
	𝑓 = range(0, stop=0.5, length=𝑓pts)
	𝑓 = ones(M)*𝑓' #2D frequency matrix for downstream algorithms
	spec = fft(sig)/(M*N/4)
	specdB = dbv.(spec[:,1:𝑓pts])

	#Convert to waveforms:
	specdB = waveform(𝑓, specdB)

	push!(plot,
		cons(:wfrm, specdB, line=set(style=:solid, color=color, width=2), label=id),
	)
	return plot
end

function plotSpec(sig; M::Int=1, id::String="signal", color=:blue)
	plot = plotSpec()
	return plotSpec!(plot, sig, id=id, color=color, M=M)
end


#==Modulator output spectrum
===============================================================================#
function plotModSpectrum()
	plot = plotSpec()
	plot.title = "Modulator Output Spectrum"
	return plot
end

function plotModSpectrum!(plot, specinfo; id::String="Simulation", color=:blue)
	𝑓 = specinfo.freq

	#Convert to waveforms:
	ePSD = waveform(𝑓, specinfo.ePSD)
	specdB = waveform(𝑓, specinfo.specdB)

	snr_str = @sprintf("SNR = %4.1fdB", specinfo.SNR) #orig. coords: (0.05,-10)
	nbw_str = @sprintf("NBW = %4.1E × 𝑓s", specinfo.NBW) #orig. coords: (0.5, -90)
	ann_str = string(snr_str, "\n", nbw_str)
	push!(plot, 
		cons(:wfrm, specdB, line=set(style=:solid, color=color, width=2), label=id),
		cons(:wfrm, ePSD, line=set(style=:solid, color=:red, width=2), label="Expected PSD"),
		cons(:atext, ann_str, y=-90, reloffset=set(x=0.5), align=:cc),
	)

	return plot
end

"""`plotModSpectrum(sig, NTF, iband, f; M=1, title, id)`

Generate plot of modulator output spectrum:
 - `sig`: Time domain ΔΣ signal.
 - `NTF`:
 - `f`: Bin location of input signal.
 - `M`: Number of modulator steps (= nlev-1).
"""
function plotModSpectrum(specinfo; id::String="Simulation", color=:blue)
	plot = plotModSpectrum()
	plotModSpectrum!(plot, specinfo, id=id, color=color)
	return plot
end


#==plotSpectrum()
===============================================================================#
function plotSpectrum()
	plot = cons(:plot, loglin, title = "Smoothed Spectrum", legend=true,
		xyaxes=set(xmin=1e-3, xmax=0.5, ymin=-120, ymax=0),
		labels=set(xaxis="Normalized Frequency", yaxis=""),
	)
	return plot
end

function plotSpectrum!(plot, X, fin; nbin::Int=8, n::Int=3, lw=1, color=:blue, id="")
	(f, p) = logsmooth(X, fin, nbin=nbin, n=n)
	pw = waveform(f, p)
	push!(plot, 
		cons(:wfrm, pw, line=set(style=:solid, color=color, width=lw), label=id),
	)
	return plot
end

"""`plotSpectrum(X, fin, nbin=8, n=3, lw=1)`

Plot a smoothed spectrum.
"""
function plotSpectrum(X, fin; nbin::Int=8, n::Int=3, lw=1, color=:blue, id="")
	plot = plotSpectrum()
	plotSpectrum!(plot, X, fin, nbin=nbin, n=n, lw=lw, color=color, id=id)
	return plot
end


#==plotExampleSpectrum
===============================================================================#
"""`plotExampleSpectrum(NTF, OSR=64, M=1, f0=0, quadrature=false, ampdB=-3, ftest=nothing, N=2^12, sizes=nothing)`

#Inputs
 - `M`: Number of modulator steps (= nlev-1).
 - `ampdB`: Test tone amplitude, relative to full-scale.
"""
function plotExampleSpectrum(NTF; OSR::Int=64, M::Int=1, f0::Float64=0, quadrature::Bool=false,
		ampdB::Float64=-3.0, ftest=nothing, N::Int=2^12, sizes=nothing)
	#Computed defaults
	if isnothing(ftest)
		ftest = default_ftest(OSR, f0=f0, quadrature=quadrature)
	end
	fband = default_fband(OSR, f0=f0, quadrature=quadrature)
	delta = 2

	usrsizes = sizes
	sizes = Dict{Symbol, Any}( #Defaults
		:lw => 1, #LineWidth
		:ms => 5, #MarkerSize
		:fs => 12, #FontSize
		:fw => "normal", #FontWeight
	)
	if isa(usrsizes, Dict) #Overwrite with user supplied values
		push!(sizes, usrsizes...)
	end

	#Plot an example spectrum
#	fin = round(Int, ftest*N)
	local u, simresult
	if !quadrature
		phi0=π/2 #Switch from sin to cos
		(u, iftest) = genTestTone_sin(ampdB, ftest, M=M, phi0=phi0, N=N)
		simresult = simulateDSM(u, NTF, nlev=M+1)
	else
		(u, iftest) = genTestTone_quad(ampdB, ftest, M=M, N=N)
		simresult = simulateQDSM(u, NTF, nlev=M+1)
	end

	plot = cons(:plot, linlin, title = "Modulator Output Spectrum", legend=true,
		labels=set(xaxis="Normalized Frequency", yaxis="PSD [dBFS/NBW]"),
		xyaxes=set(ymin=-140, ymax=0),
	)

	specinfo = calcSpecInfo(simresult.v, NTF, fband, ftest, M=M, quadrature=quadrature)
	if !quadrature
		if 0==f0
#			set(plot, loglin, xyaxes=set(xmin=.001))
		end

		#circular smoothing; not re-scaled
		spec_smoothed = circ_smooth(abs.(specinfo.spec0).^2, 16)
		A = dbv(specinfo.spec0[specinfo.itest+1])

		#Convert to waveforms & add to plot:
		𝑓 = specinfo.freq; 𝑓pts = length(𝑓)
		specw = waveform(𝑓, specinfo.specdB)
		spec_sw = waveform(𝑓, dbp.(spec_smoothed[1:𝑓pts]))
		ePSDw = waveform(𝑓, specinfo.ePSD)

		push!(plot, 
			cons(:wfrm, specw, line=set(style=:solid, color=:cyan, width=2), label=""),
			cons(:wfrm, spec_sw, line=set(style=:solid, color=:blue, width=2), label=""),
			cons(:wfrm, ePSDw, line=set(style=:solid, color=:magenta, width=1), label=""),
		)
		alignleft = f0<0.25
		align = alignleft ? :tl : :tr
		annx  = alignleft ? f0+0.5/OSR : f0-0.5/OSR
		tstr = @sprintf(" SQNR = %.1fdB\n @ A = %.1fdBFS\n  & OSR = %.0f\n",
			specinfo.SNR, A, OSR
		)
		nbwstr = @sprintf("NBW=%.1e ", specinfo.NBW)
		#TODO: use `sizes` atext:
		push!(plot, 
			cons(:atext, tstr, x=annx, y=-5, align=align),
			cons(:atext, nbwstr, x=0.5, y=-135, align=:cr),
		)
	else
		error("plotExampleSpectrum not implemented for quadrature modulators.")
#=
		spec0 = fftshift(spec0/2)
		freq = range(-0.5, stop=0.5, length=N+1)
		freq = deleteat!(collect(freq), N+1)
		plot(freq,dbv(spec0),'c','Linewidth',1);
		spec_smoothed = circ_smooth(abs(spec0).^2, 16);
		plot(freq, dbp(spec_smoothed), 'b', 'Linewidth', 3);
		Snn = abs(evalTF(NTF,exp(2i*pi*freq))).^2 * 2/12 * (delta/M)^2;
		plot(freq, dbp(Snn*NBW), 'm', 'Linewidth', 1);
		SNR = calculateSNR(spec0(N/2+1+[f1_bin:f2_bin]),fin-f1_bin);
		msg = sprintf('SQNR = %.1fdB\n @ A=%.1fdBFS & OSR=%.0f\n', ...
		SNR, dbv(spec0(N/2+fin+1)), OSR );
		if f0>=0
			text(f0-0.05, -15, msg,'Hor','Right');
		else
			text(f0+0.05, -15, msg,'Hor','Left');
		end
		text(-0.5,-135,sprintf(' NBW=%.1e',NBW),'hor','left');
		figureMagic([-0.5 0.5],0.125,2, [-140 0],10,2);
=#
	end

	return plot
end


"""`plotExampleSpectrum(dsm, NTF=nothing; ampdB=-3, ftest=nothing, N=2^12, sizes=sizes)`
"""
function plotExampleSpectrum(dsm::AbstractDSM, NTF=nothing; ampdB=-3.0, ftest=nothing, N=2^12, sizes=nothing)
	isnothing(NTF) && (NTF = synthesizeNTF(dsm))
	return plotExampleSpectrum(NTF, OSR=dsm.OSR, M=dsm.M, f0=dsm.f0, quadrature=isquadrature(dsm),
		ampdB=ampdB, ftest=ftest, N=N, sizes=sizes
	)
end

#Last line
