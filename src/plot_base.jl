#RSDeltaSigmaPort: Base plot generation tools
#-------------------------------------------------------------------------------


#==Useful constants
===============================================================================#
const linlin = cons(:a, xyaxes=set(xscale=:lin, yscale=:lin))
const linlog = cons(:a, xyaxes=set(xscale=:lin, yscale=:log))
const loglin = cons(:a, xyaxes=set(xscale=:log, yscale=:lin))
const loglog = cons(:a, xyaxes=set(xscale=:log, yscale=:log))

const dfltline = cons(:a, line=set(style=:solid, color=:blue, width=2))


#==Waveform builders
===============================================================================#

#-------------------------------------------------------------------------------
function _wfrm__(x::Vector, y::Vector, sweepid::String)
	validatexy_vec(x, y)
	return DataF1(x,y)
end
function _wfrm__(x::Array{Float64,2}, y::Array{Float64,2}, sweepid::String)
	validatexy_1param(x, y)
	N = size(x, 1)
	wfrm = fill(DataRS{DataF1}, PSweep(sweepid, collect(1:N))) do i
		return DataF1(x[i,:], y[i,:]) #Returns from implicit "do" function
	end
	return wfrm
end

"""`waveform(x, y; sweepid="i")`

Create a waveform object with a given string to identify the sweep.
"""
waveform(x, y; sweepid::String="i") = _wfrm__(x, y, sweepid)

#-------------------------------------------------------------------------------
function _stairs__(x::Vector, y::Vector, sweepid::String)
	validatexy_vec(x, y)
	if length(x) < 1
		throw("Does not support zero-length arrays")
	end
	xs = similar(x, 2*length(x))
	ys = similar(y, 2*length(y))

	ii = io = 1
	while ii < length(x)
		xs[io] = x[ii]
		xs[io+1] = x[ii+1]
		ii += 1; io +=2
	end
	if ii <= length(x) #last point
		Δx = 1 #Step size (if x has a single point)
		if length(x) > 1; Δx = x[end]-x[end-1]; end
		xs[io] = x[ii]
		xs[io+1] = x[ii]+Δx #Use final step size
	end

	ii = io = 1
	while ii <= length(y)
		_y = y[ii]
		ys[io] = _y
		ys[io+1] = _y
		ii += 1; io +=2
	end
	return DataF1(xs, ys)
end

function _stairs__(x::Array{Float64,2}, y::Array{Float64,2}, sweepid::String)
	validatexy_1param(x, y)
	N = size(x, 1)
	wfrm = fill(DataRS{DataF1}, PSweep(sweepid, collect(1:N))) do i
		return _stairs__(x[i,:], y[i,:], "") #Returns from implicit "do" function
	end
	return wfrm
end

"""`wfrm_stairs(x, y; sweepid="i")`

Create a staircase-waveform object with a given string to identify the sweep.
"""
wfrm_stairs(x, y; sweepid::String="i") = _stairs__(x, y, sweepid)


#==Misc. draw functions
===============================================================================#
function plot_unitcircle!(plot; color=:black, width=2, nseg::Int=360)
	Θ = range(0, stop=2π, length=nseg+1)
	cplxcircle = exp.((2π*j)*Θ)
	circle = DataF1(real(cplxcircle), imag(cplxcircle))
	push!(plot,
		cons(:wfrm, circle, line=set(style=:solid, color=color, width=width))
	)
	return plot
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


#==Pole-zero plots
===============================================================================#
#Create empty pole-zero plot
function plotPZ()
	plot = cons(:plot, title = "Pole-Zero Plot", legend=false,
		xyaxes=set(xscale=:lin, yscale=:lin, xmin=-1.1, xmax=1.1, ymin=-1.1, ymax=1.1),
		labels=set(xaxis="Real", yaxis="Imag"),
	)
	plot_unitcircle!(plot)
	return plot
end

"""`function plotPZ!(H,color=:black,markersize=5,list=false)`
Plot the poles and zeros of a transfer function.
If list=true, a list of the poles and zeros is superimposed on the plot.
"""
function plotPZ!(plot, H::ZPKData; color=:blue, markersize=2, list::Bool=false)
	isVector = isa(color, Vector)
	color_pole = isVector ? color[1] : color
	color_zero = isVector ? color[end] : color
	attrib_poles = cons(:a,
		glyph=set(shape=:x, color=color_pole, size=markersize), line=set(style=:none, width=2)
	)
	attrib_zeros = cons(:a,
		glyph=set(shape=:o, color=color_zero, size=markersize), line=set(style=:none, width=2)
	)
	
	z, p, k = _zpkdata(H)

	# Plot x and o for poles and zeros, respectively:
	push!(plot, cons(:wfrm, DataF1(real(p), imag(p)), attrib_poles))
	if !isempty(z) #Add zeros
		push!(plot, cons(:wfrm, DataF1(real(z), imag(z)), attrib_zeros))
	end

#=
	if list
		# List the poles and zeros
		p = cplxpair(p);
		y = 0.05*(ceil(length(p)/2)+1);
		str = 'Poles:               ';
		text( 0, y, str, 'Hor', 'Right', 'Ver', 'Mid'); y = y - 0.1;
		for i = 1:2:length(p);
			if abs(imag(p(i))) < 1e-6
				str = sprintf('%+.4f      ', real(p(i)) );
			else
				str = sprintf( '%+.4f+/-j%.4f  ', real(p(i)), abs(imag(p(i))) );
			end
			text( 0, y, str, 'Hor', 'Right', 'Ver', 'Mid'); y = y - 0.1;
		end
		if !isempty(z)
			z = z( !isnan(z) && !isinf(z) );
			z = cplxpair(z);
			y = 0.05*(ceil(length(z)/2)+1);
			str = '        Zeros:';
			text( 0, y, str, 'Hor', 'Left', 'Ver', 'Mid'); y = y - 0.1;
			for i = 1:2:length(z);
				if abs(imag(z(i))) < 1e-6
					str = sprintf('%+.4f      ', real(z(i)) );
				else
					str = sprintf( '  %+.4f+/-j%.4f', real(z(i)), abs(imag(z(i))) );
				end
				text( 0, y, str, 'Hor', 'Left', 'Ver', 'Mid');
				y = y - 0.1;
			end
		end
	end
=#
	return plot
end

plotPZ(H::ZPKData, args...; kwargs...) = plotPZ!(plotPZ(), H::ZPKData, args...; kwargs...)


#==NTF plots
===============================================================================#
function plotNTF(isbandpass)
	dBmin = -100 #Limit dB window to finite value
	𝑓min_log = 10^-2 #Limit min 𝑓 for log(𝑓-norm) plots
	nogrid = set(vmajor = false, vminor = false, hmajor=false, hminor=false)
	vgridonly = set(vmajor = true, vminor = true, hmajor=false, hminor=false)

	plot_PZ = plotPZ()
	plot_PZ.title = "NTF Poles and Zeros"

	plot_NTFlin = cons(:plot, legend=false,
		title = "NTF Magnitude Response",
		xyaxes=set(xscale=:lin, yscale=:lin, ymin=dBmin),
		labels=set(xaxis="Normalized frequency (1→f_s)", yaxis="dB"),
	)

	plot_NTFlog = cons(:plot, legend=false, grid = vgridonly,
#		title = "NTF Magnitude Response",
		xyaxes=set(xscale=:log, yscale=:lin, xmin=𝑓min_log, ymin=dBmin),
		labels=set(xaxis="Normalized frequency (1→f_B)", yaxis="dB"),
	)
	plot_NTFzoom = cons(:plot, legend=false, grid = vgridonly,
#		title = "NTF Magnitude Response",
		xyaxes=set(xscale=:lin, yscale=:lin, xmin=-0.6, ymax=0.6, ymin=dBmin),
		labels=set(xaxis="Normalized frequency offset", yaxis="dB"),
	)
	if !isbandpass; plot_NTFzoom = plot_NTFlog; end

	pcoll = push!(cons(:plot_collection, title="Sample Plot"),
		plot_PZ, plot_NTFlin, plot_NTFzoom
	)

	pcoll.bblist = [
		BoundingBox(0, 0.5, 0, 1), #Pole-zero
		BoundingBox(0.5, 1, 0, 0.5), #NTF linf
		BoundingBox(0.5, 1, 0.5, 1), #NTF zoom
	]

	return pcoll
end

function plotNTF!(pcoll, NTF::ZPKData, OSR::Int; color=:blue, f0::Real=0, STF=nothing)
	dfltline = cons(:a, line=set(style=:solid, color=color, width=2))
	gainline = cons(:a, line=set(style=:solid, color=:red, width=2))
	markerline = cons(:a, line=set(style=:dash, width=2.5))
	plot_PZ, plot_NTFlin, plot_NTFzoom = pcoll.plotlist
	isbandpass = f0!=0

	#Plot poles and zeros
	plotPZ!(plot_PZ, NTF, color=color)

	#Plot full-range NTF
	𝑓 = vcat(range(0, stop=0.75/OSR, length=100), range(0.75/OSR, stop=0.5, length=100))
	if isbandpass
		f1 = f0-1/(2*OSR); f2 = f0+1/(2*OSR)
		𝑓 = vcat(range(0, stop=f1, length=50), range(f1, stop=f2, length=100), range(f2, stop=0.5, length=50))
	end
	z = exp.(2j*π*𝑓)
	magNTF = dbv(evalTF(NTF, z))
	push!(plot_NTFlin, 
		cons(:wfrm, DataF1(𝑓, magNTF), dfltline, label="|NTF|")
	)
	if STF != nothing
		plot_NTFlin.title = "NTF/STF Magnitude Response"
		magSTF = dbv(evalTF(STF,z))
		push!(plot_NTFlin, 
			cons(:wfrm, DataF1(𝑓, magSTF), gainline, label="|STF|")
		)
	end


	#Plot "zoomed-in" NTF
	𝑓start = 0.01
	𝑓 = collect(range(𝑓start, stop=1.2, length=200)/(2*OSR))
	∫rng = [0, 0.5/OSR]
	if isbandpass
		f1 = f0-0.3/OSR; f2 = f0+0.3/OSR
		𝑓 = collect(range(f1, stop=f2, length=100))
		∫rng = f0 .+ [-.25, .25]/OSR
	end
	z = exp.(2j*π*𝑓)
	𝑓norm = (𝑓 .- f0)*(2*OSR)
	magNTF = dbv(evalTF(NTF, z))
	σ_NTF = dbv(rmsGain(NTF, ∫rng[1], ∫rng[2]))
	rmsgain_str = @sprintf("RMS Gain = %5.0fdB", σ_NTF)
	push!(plot_NTFzoom, 
		cons(:wfrm, DataF1(𝑓norm, magNTF), dfltline, label="|NTF|"),
		cons(:atext, rmsgain_str, y=σ_NTF, offset=set(y=3), reloffset=set(x=0.5), align=:bc),
		cons(:hmarker, σ_NTF, markerline, strip=2),
	)

	return pcoll
end

plotNTF(NTF::ZPKData, args...; f0=0, kwargs...) =
	plotNTF!(plotNTF(f0!=0), NTF::ZPKData, args...; f0=f0, kwargs...)


#==Time domain plots of modulator
===============================================================================#
function plotModTransient(inputSig, outputSig, otherSig...; legend::Bool=true, color=:blue)
	plot = cons(:plot, linlin, title = "Modulator Input & Output", legend=legend,
		labels=set(xaxis="Sample Number", yaxis="Amplitude [V]"),
	)

	#Ensure we have 2D-Array{}s, not Vector{}s:
	inputSig = conv2seriesmatrix2D(inputSig)
	outputSig = conv2seriesmatrix2D(outputSig)

	#Generate matrix of sample numbers:
	M = size(inputSig,1) #Number of rows (ex: quantizers)
	N = size(inputSig,2) #Number of samples
	n = 1:N
	sn = ones(Float64, M)*n' #Matrix of "sample number" (Use Float64s for downstream functions)

	inputSig = wfrm_stairs(sn, inputSig)
	outputSig = wfrm_stairs(sn, outputSig)

	push!(plot,
		cons(:wfrm, outputSig, line=set(style=:solid, color=color, width=2), label="output"),
		cons(:wfrm, inputSig, line=set(style=:solid, color=:red, width=2), label="input"),
#		cons(:wfrm, otherSig, line=set(style=:solid, color=:black, width=2), label="y"),
	)
	return plot
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


#==SNR plot
===============================================================================#
function plotSNR(; legend::Bool=true,
		title::String="SNR: Theory vs Simulation"
	)

	plot = cons(:plot, linlin, title=title, legend=legend,
		#xyaxes=set(xmin=-100, xmax=0, ymin=0, ymax=100),
		xyaxes=set(ymax=100), #Clip infinite gains
		labels=set(xaxis="Input Level (dBFS)", yaxis="SQNR (dB)"),
	)

	return plot
end

function plotSNR!(plot, sig, NTF, OSR::Int; f0::Float64=0.0, nlev::Int=2,
		id::String="simulation", color=:green
	)
	sig = conv2seriesmatrix2D(sig)
	M, N = size(sig)

	snr, amp = simulateSNR(NTF, OSR, f0=f0, nlev=nlev)
	pk_snr, pk_amp = peakSNR(snr, amp) #Single values

	snr_pred, amp_pred = nothing, nothing
	if nlev == 2 #predictSNR() only supports 2 levels
		snr_pred, amp_pred = predictSNR(NTF, OSR, f0=f0)
		amp_pred = conv2seriesmatrix2D(amp_pred*1.0) #Ensure amplitudes are Float64
	end

	#Convert to waveforms:
	snr = waveform(amp, snr)

	simglyph = cons(:a, glyph=set(shape=:o, size=1, color=color, fillcolor=color))

	snr_str = @sprintf("peak SNR = %4.1fdB\n@ OSR = %d", pk_snr, OSR) #orig. coords: (-25, 85)
	push!(plot, 
		cons(:wfrm, snr, simglyph, line=set(style=:dashdot, color=color, width=2), label="simulation"),
		cons(:atext, snr_str, x=-25, y=85, align=:cr),
	)
	if !isnothing(snr_pred)
		snr_pred = waveform(amp_pred, snr_pred) #Convert to waveform
		push!(plot,
			cons(:wfrm, snr_pred, line=set(style=:solid, color=:red, width=2), label="theory"),
		)
	end
	return plot
end

function plotSNR(sig, NTF, OSR::Int; f0::Float64=0.0, nlev::Int=2, legend::Bool=true,
		title::String="SNR: Theory vs Simulation", id::String="simulation", color=:green
	)
	plot = plotSNR(legend=legend, title=title)
	plotSNR!(plot, sig, NTF, OSR, f0=f0, nlev=nlev, id=id, color=color)
	return plot
end


#==Plot state maxima
===============================================================================#
function plotStateMaxima(u, ABCD; nlev::Int=2, N::Int=10^5)
	title = "Simulated State Maxima"; color = :blue
	T = ones(1,N)
	nu, nq, order = get_nu_nq_order(T, ABCD, nlev)
	u = collect(u)
	maxima = zeros(order,length(u))

	for i = 1:length(u)
		ui = u[i]
		v,xn,xmax = simulateDSM(ui*T, ABCD, trackmax=true)
		maxima[:,i] = xmax[:]
		if any(xmax .> 1e2)
			umax = ui
			u = u[1:i]
			maxima = maxima[:,1:i]
			break
		end
	end

	plot = cons(:plot, linlog, title = title, legend=true,
#		xyaxes=set(xmin=0, xmax=0.6, ymin=1e-4, ymax=10),
		xyaxes=set(ymin=1e-4, ymax=10),
		labels=set(xaxis="DC input", yaxis="Peak value"),
	)
	for i in 1:order
		id = "state $i"
		_maxima = waveform(u, maxima[i,:])
		push!(plot,
			cons(:wfrm, _maxima, line=set(style=:solid, width=2), label=id),
		)
	end
	return plot
end


#Last line
