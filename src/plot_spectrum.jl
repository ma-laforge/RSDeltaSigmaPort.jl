#RSDeltaSigmaPort: Generate frequency spectrum plots
#-------------------------------------------------------------------------------

#==Helper functions
===============================================================================#
"""`(f1, f2) = ds_f1f2(OSR=64; f0=0, iscomplex=0)`
"""
function ds_f1f2(OSR::Int; f0::Float64=0, iscomplex=false)
	local f1, f2
	if iscomplex
		f1 = f0-0.5/OSR
		f2 = f0+0.5/OSR
	else
		if f0>0.25/OSR
			f1 = f0-0.25/OSR
			f2 = f0+0.25/OSR
		else
			f1 = 0
			f2 = 0.5/OSR
		end
	end
	return (f1, f2)
end


#==Frequency-spectrum plots
===============================================================================#
function plotSpec()
	plot = cons(:plot, linlin, title = "Output Spectrum", legend=true,
		xyaxes=set(xmin=0, xmax=0.5, ymin=-120, ymax=0),
		labels=set(xaxis="Normalized Frequency", yaxis="dBFS"),
	)
	return plot
end

function plotSpec!(plot, sig; id::String="signal", color=:blue, sp2p::Float64=1.0)
	sig = conv2seriesmatrix2D(sig)
	M, N = size(sig)
	Nspec = div(N,2)+1 #Only need half the spectrum
	𝑓 = range(0, stop=0.5, length=Nspec)
	𝑓 = ones(M)*𝑓' #2D frequency matrix for downstream algorithms
	spec = fft(sig)/(sp2p*N/4)
	specdB = dbv.(spec[:,1:Nspec])

	#Convert to waveforms:
	specdB = waveform(𝑓, specdB)

	push!(plot,
		cons(:wfrm, specdB, line=set(style=:solid, color=color, width=2), label=id),
	)
	return plot
end

function plotSpec(sig; id::String="signal", color=:blue, sp2p::Float64=1.0)
	plot = plotSpec()
	return plotSpec!(plot, sig, id=id, color=color, sp2p=sp2p)
end


#==Modulator output spectrum
===============================================================================#
function plotModSpectrum(;title::String="Modulator Output Spectrum")
	plot = plotSpec()
	plot.title = title
	return plot
end

function plotModSpectrum!(plot, sig, NTF, iband::IndexRange, f::Int; sp2p::Float64=1.0,
		id::String="simulation", color=:blue
	)
	sig = conv2seriesmatrix2D(sig)

	M, N = size(sig)
	if M!=1
		throw("FIXME: replicate hann window rows for M>0")
	end

	#Compute windowed signal spectrum & SNR
	swnd = sig .* ds_hann(N)' #windowed (1xN array)
	spec = fft(swnd)/(sp2p*N/4)
	snr = calculateSNR(spec[iband], f) #Must integrate from entire spectrum.

	#Compute expected PSD
	NBW = 1.5/N
	Nspec = div(N,2)+1 #Only need half the spectrum for display purposes
	𝑓 = collect(range(0, stop=0.5, length=Nspec))
	Sqq = 4 * evalTF(NTF,exp.(2j*pi*𝑓)).^2 / (3*sp2p^2)
	psde = dbp.(Sqq*NBW)

	#Convert to waveforms:
	psde = waveform(𝑓, psde)

	plot = plotSpec!(plot, swnd, id=id, color=color, sp2p=sp2p)

	snr_str = @sprintf("SNR = %4.1fdB", snr) #orig. coords: (0.05,-10)
	nbw_str = @sprintf("NBW = %4.1E × 𝑓s", NBW) #orig. coords: (0.5, -90)
	ann_str = string(snr_str, "\n", nbw_str)
	push!(plot, 
		cons(:wfrm, psde, line=set(style=:solid, color=:red, width=2), label="Expected PSD"),
		cons(:atext, ann_str, y=-90, reloffset=set(x=0.5), align=:cc),
	)

	return plot
end

"""`plotModSpectrum(sig, NTF, iband, f; sp2p=1.0, title, id)`

Generate plot of modulator output spectrum:
 - sig: Time domain ΔΣ signal.
 - NTF:
 - f: Bin location of input signal.
 - sp2p: Peak-to-peak amplitude (V) of input signal.
"""
function plotModSpectrum(sig, NTF, iband::IndexRange, f::Int; sp2p::Float64=1.0,
		title::String="Modulator Output Spectrum", id::String="simulation", color=:blue
	)
	plot = plotModSpectrum(title=title)
	plotModSpectrum!(plot, sig, NTF, iband, f, sp2p=sp2p, id=id, color=color)
	return plot
end


#==plotExampleSpectrum
===============================================================================#
"""`plotExampleSpectrum(ntf, osr=64, M=1, f0=0, quadrature=false, Atest=-3, Ftest=nothing, N=2^12, sizes=nothing)`
"""
function plotExampleSpectrum(ntf; osr=64, M=1, f0=0, quadrature=false,
		Atest=-3, Ftest=nothing, N=2^12, sizes=nothing)

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

	f1, f2 = ds_f1f2(osr, f0=f0, iscomplex=quadrature)

	#Computed defaults
	if nothing == Ftest
		if 0==f0 || quadrature
			Ftest = 0.15/osr
		else
			Ftest = f0 + 0.08/osr
		end
	end

	delta = 2

	#Plot an example spectrum
	Amp = undbv(Atest) #Test tone amplitude, relative to full-scale.
	f1_bin = round(Int, f1*N)
	f2_bin = round(Int, f2*N)
	fin = round(Int, Ftest*N)
	local u, v
	if !quadrature
		u = Amp*M*cos.((2π/N)*fin * (0:N-1))
		(v,) = simulateDSM(u, ntf, nlev=M+1)
	else
		u = Amp*M*exp.((2π*j/N)*fin * (0:N-1))
		(v,) = simulateQDSM(u, ntf, nlev=M+1)
	end

	plot = cons(:plot, linlin, title = "Modulator Output Spectrum", legend=true,
		labels=set(xaxis="Normalized Frequency", yaxis="PSD [dBFS/NBW]"),
		xyaxes=set(ymin=-140, ymax=0),
	)

	NBW = 1.5/N
	sp2p = M #peak-to-peak signal swing
	sigwnd = v .* ds_hann(N)' #windowed signal (1xN array)
	spec0 = fft(sigwnd)/(sp2p*N/4)
	if !quadrature
		fpts = div(N,2)+1
		freq = range(0, stop=0.5, length=fpts)
		if 0==f0
			set(plot, loglin, xyaxes=set(xmin=.001))
		end
		#circular smoothing; not re-scaled
		spec_smoothed = circ_smooth(abs.(spec0).^2, 16)
		Snn = abs.(evalTF(ntf, exp.(2π*j*freq))).^2 * 2*2/12 * (delta/M)^2
		snr = calculateSNR(spec0[f1_bin+1:f2_bin+1], fin-f1_bin)

		#Convert to waveforms & add to plot:
		specw = waveform(freq, dbv.(spec0[1:fpts]))
		spec_sw = waveform(freq, dbp.(spec_smoothed[1:fpts]))
		Snnw = waveform(freq, dbp.(Snn*NBW))
		push!(plot, 
			cons(:wfrm, specw, line=set(style=:solid, color=:cyan, width=2), label=""),
			cons(:wfrm, spec_sw, line=set(style=:solid, color=:blue, width=2), label=""),
			cons(:wfrm, Snnw, line=set(style=:solid, color=:magenta, width=1), label=""),
		)
		alignleft = f0<0.25
		align = alignleft ? :tl : :tr
		annx  = alignleft ? f0+0.5/osr : f0-0.5/osr
		tstr = @sprintf(" SQNR = %.1fdB\n @ A = %.1fdBFS\n  & OSR = %.0f\n",
			snr, dbv(spec0[fin+1]), osr
		)
		nbwstr = @sprintf("NBW=%.1e ", NBW)
		#TODO: use `sizes` atext:
		push!(plot, 
			cons(:atext, tstr, x=annx, y=-5, align=align),
			cons(:atext, nbwstr, x=0.5, y=-135, align=:cr),
		)
	else
#=
		spec0 = fftshift(spec0/2)
		freq = range(-0.5, stop=0.5, length=N+1)
		freq = deleteat!(collect(freq), N+1)
		plot(freq,dbv(spec0),'c','Linewidth',1);
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
=#
	end

	return plot
end


"""`plotExampleSpectrum(dsm, ntf=nothing; sizes=sizes, Atest=-3, Ftest=nothing, N=2^12)`
"""
function plotExampleSpectrum(dsm::AbstractDSM, ntf=nothing; sizes=nothing, Atest=-3, Ftest=nothing, N=2^12)
	isnothing(ntf) && (ntf = synthesizeNTF(dsm))
	return plotExampleSpectrum(ntf, osr=dsm.osr, M=dsm.M, f0=dsm.f0, quadrature=isquadrature(dsm),
		Atest=Atest, Ftest=Ftest, N=N, sizes=sizes
	)
end

#Last line
