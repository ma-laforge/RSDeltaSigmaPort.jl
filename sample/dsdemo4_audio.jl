#dsdemo4_audio.jl: A GUI-based demonstration of delta-sigma with sound (Currently GUI-less).
using RSDeltaSigmaPort

@warn("dsdemo4_audio requires additional libraries:\n - MAT.jl\n - WAV.jl")

if !@isdefined(DSAudioDemo) #Hack to avoid continuously appending docs
@doc """# `DSAudioDemo`: Delta-Sigma Audio Demo Module

 - Set source, modulator and decimation filter paramters,
 - then click "Go."
 - Click on input or output waveforms to listen to them.
 - Click on the output spectrum to listen to the amplified error.

## TODO
 - Implement GUI, as originally intended.
""" DSAudioDemo
end

module DSAudioDemo

using RSDeltaSigmaPort
using RSDeltaSigmaPort.EasyPlot #set, cons
using RSDeltaSigmaPort: fft
using RSDeltaSigmaPort: linlin, loglin
import Printf: @sprintf
import MAT, WAV


#==Constants
===============================================================================#
const 𝑓ₛout = 8192 #Output sample rate is fixed
const j=im

"ABCD matrices of pre-defined modulators:"
const MODABCD_PRESETS = Dict{Int, Array}(
	1 => [1 1 -1; 1 0 0],
	2 => [1 0 1 -1; 1 1 1 -2; 0 1 0 0],
)


#==Types
===============================================================================#
abstract type AbstractSource; end

struct SrcRamp <: AbstractSource; end

mutable struct SrcRaw <: AbstractSource
	u::Vector #{Float64} #Oversampled signal
	u0::Vector #{Float64} #Decimated signal
end

mutable struct SrcSine <: AbstractSource
	f::Int #Frequency (Hz)
	A::Float64
end
SrcSine(;f=500, A=0.5) = SrcSine(Int(f), A)

struct Modulator
	sys
	fsout::Int #Output sample rate
	OSR::Int
	sincorder::Int
end


#==Constructors
===============================================================================#
function _Modulator(sys, 𝑓ₒₛ::Int, OSR::Nothing, sincorder::Int)
	#Round 𝑓ₒₛ to a multiple of 𝑓ₛout:
	OSR = round(Int, 𝑓ₒₛ/𝑓ₛout)
	OSR = max(1, OSR)
	𝑓ₒₛnew = 𝑓ₛout*OSR #𝑓ₛout is always hardcoded in dsdemo4

	if 𝑓ₒₛ < 𝑓ₛout
		@warn("Requested (over)-sampling frequency, 𝑓ₒₛ < 𝑓ₛout." *
			"\nSetting 𝑓ₒₛ = 𝑓ₛout = $𝑓ₛout Hz."
		)
	elseif 𝑓ₒₛ != 𝑓ₒₛnew
		@warn("Adjusting (over)-sampling frequency\nfrom 𝑓ₒₛ = $𝑓ₒₛ" *
			"\nto 𝑓ₒₛ = $𝑓ₒₛnew (OSR = $OSR)."
		)
	end

	return Modulator(sys, 𝑓ₛout, OSR, sincorder)
end
_Modulator(sys, 𝑓ₒₛ::Nothing, OSR::Int, sincorder::Int) =
	Modulator(sys, 𝑓ₛout, OSR, sincorder) #𝑓ₛout is always hardcoded in dsdemo4
_Modulator(sys, 𝑓ₒₛ::Nothing, OSR::Nothing, sincorder::Int) =
	_Modulator(sys, 𝑓ₒₛ, 32, sincorder) #Default OSR
_Modulator(sys, 𝑓ₒₛ, OSR, sincorder::Int) =
	throw(ArgumentError("Modulator: Invalid type for 𝑓ₒₛ or OSR (or overspecified)."))

Modulator(preset; fos=nothing, OSR=nothing, sincorder::Int=2) =
	_Modulator(MODABCD_PRESETS[preset], fos, OSR, sincorder)


#==Helper functions/ data generators
===============================================================================#
function validate(m::Modulator)
	(m.OSR < 1) &&
		throw("Modulator: Invalid OSR value: $(m.OSR).")
	return
end

function generate(src::SrcSine, m::Modulator; N=10)
	(src.f > m.fsout/2) &&
		throw("Source 𝑓 = $(src.f) too high. Must be <= 𝑓ₛout/2.")
	𝑓ₒₛ = m.fsout*m.OSR
	𝜔₀ = 2π*src.f
	tidx = (0:N-1); t = tidx/𝑓ₒₛ
	u = src.A * sin.((𝜔₀/𝑓ₒₛ) * tidx) .* ds_hann(N)
	u0 = u[1:m.OSR:end]
	return (t, u, u0)
end

function generate(src::SrcRamp, m::Modulator; N=10)
	𝑓ₒₛ = m.fsout*m.OSR
	tidx = (0:N-1); t = tidx/𝑓ₒₛ
	u = collect(range(-0.7, stop=0.7, length=N))
	u0 = u[1:m.OSR:end]
	return (t, u, u0)
end

function generate(src::SrcRaw, m::Modulator; N=10)
	𝑓ₒₛ = m.fsout*m.OSR
	N = length(src.u)
	tidx = (0:N-1); t = tidx/𝑓ₒₛ
	return (t, src.u, src.u0)
end

function play(data::Vector)
	WAV.wavplay(data, 𝑓ₛout)
	return :PLAYBACK_COMLETE
end


#==Main algorithms
===============================================================================#
function run(m::Modulator; Tsim::Float64=2.0, input=SrcSine(f=500))
	validate(m)
	(m.fsout != 𝑓ₛout) && throw("dsdemo4 ony supports modulation schemes where 𝑓ₛout = $𝑓ₛout.")
	fos = 𝑓ₛout*m.OSR #(over)-sampling frequency

	#Compute signal `u` depending on desired source:
	N = round(Int, Tsim*fos) #Unless overwritten
	(t, u, u0) = generate(input, m, N=round(Int, Tsim*fos))
	OSR = m.OSR #Copy for string interpolation

	#Plot signal `u0`
	plot_tdec = cons(:plot, linlin, title="Decimated Signal (dec=$OSR)", legend=false,
		xyaxes=set(ymin=-1, ymax=1),
		labels=set(xaxis="Time [s]", yaxis="u₀(t), w(t)"),
	)

	N = length(u0)
	t = (0:N-1)/𝑓ₛout
	u0w = waveform(t, u0)
	push!(plot_tdec,
		cons(:wfrm, u0w, line=set(style=:solid, color=:blue, width=2), label="u₀"),
	)

	#Plot U(𝑓), from 0 to 𝑓ₛout/2
	plot_𝑓dec = cons(:plot, linlin, title="Spectrum (dec=$OSR)", legend=false,
		xyaxes=set(xmin=0, xmax=𝑓ₛout/2, ymin=-160, ymax=0),
		labels=set(xaxis="Frequency [Hz]", yaxis="U₀(𝑓), W(𝑓)"),
	)

	if typeof(input) in [SrcSine, SrcRaw]
		N = length(u0)
		local U
		if isa(input, SrcSine)
			U = fft(u0)/(N/4)
		else #SrcRaw
			U = fft(applywnd(u0, ds_hann(N)))/(N/4)
		end
		𝑓 = range(0, stop=𝑓ₛout, length=N+1); 𝑓 = 𝑓[1:div(N,2)+1]
		UdB = dbv.(U[keys(𝑓)])
		UdBw = waveform(𝑓, UdB)
		push!(plot_𝑓dec,
			cons(:wfrm, UdBw, line=set(style=:solid, color=:blue, width=2), label="U₀"),
		)
	end

	#Plot ΔΣ signals (time domain):
	Nplot = 300 #Number of samples to plot
	simresult = simulateDSM(u, m.sys)
	v = simresult.v; y = simresult.y
	q = v - y #Quantization error; assumes quantizer gain = 1.
	N = length(v)
	n = 1:N
	if N>Nplot #Sample values from the middle:
		n = floor.(Int, N/2-Nplot/2:N/2+Nplot/2-1)
	end

	plot_tos = cons(:plot, linlin, title="Signals @ Modulator Input/Output", legend=false,
		xyaxes=set(ymin=-1.2, ymax=1.2),
		labels=set(xaxis="samples", yaxis="v(t)"),
	)

	_x = collect(0:length(n)-1)
	vw = wfrm_stairs(_x, v[n])
	uw = waveform(_x, u[n])

	push!(plot_tos,
		cons(:wfrm, vw, line=set(style=:solid, color=:blue, width=2), label="v"),
		cons(:wfrm, uw, line=set(style=:solid, color=:green, width=2), label="u"),
	)

	#Plot V(𝑓), from 0 to 𝑓ₒₛ/2. Use the middle Nfft points of v
	N = length(v)
	Nfft = min(N, 16*8192)
	n = array_round((N-Nfft)/2+1):array_round((N+Nfft)/2)
	V = fft(applywnd(v[n], ds_hann(Nfft)))/(Nfft/4)
	inBin = ceil(Int, Nfft/1000)
	if isa(input, SrcSine)
		inBin = round(Int, input.f/fos*Nfft+1) #Bin of tone
	end
	(𝑓nrm, Vp) = logsmooth(V, inBin)

	plot_𝑓os = cons(:plot, loglin, title="Spectrum (OSR=$OSR)", legend=false,
		xyaxes=set(xmin=100, xmax=fos/2, ymin=-160, ymax=0),
		labels=set(xaxis="Frequency [Hz]", yaxis="V(𝑓)"),
	)

	Vpw = waveform(𝑓nrm*fos, Vp)
	nbw_str = @sprintf("NBW = %.1f Hz", fos*1.5/Nfft) #orig. coords: (Fs/2, -90); :cr
	push!(plot_𝑓os,
		cons(:wfrm, Vpw, line=set(style=:solid, color=:blue, width=2), label="V"),
		cons(:atext, nbw_str, y=-90, reloffset=set(x=0.95), align=:cr),
	)

	#Compute w
	w = sinc_decimate(v, m.sincorder, m.OSR)
	filtered_q = sinc_decimate(q, m.sincorder, m.OSR)
	N = length(w)
	t = collect(0:N-1)/𝑓ₛout
	ww = wfrm_stairs(t, w)
	push!(plot_tdec,
		cons(:wfrm, ww, line=set(style=:solid, color=:red, width=2), label="w"),
	)

	#Plot W(𝑓), from 0 to 𝑓ₛout/2
	if typeof(input) in [SrcSine, SrcRaw]
		Nfft = length(w)
		local W
		if isa(input, SrcSine)
			W = fft(w)/(N/4)
		else
			W = fft(applywnd(w, ds_hann(N)))/(N/4)
		end
		𝑓 = range(0, stop=𝑓ₛout, length=Nfft+1); 𝑓 = 𝑓[1:div(Nfft,2)+1]
		WdB = dbv.(W[keys(𝑓)])
		WdBw = waveform(𝑓, WdB)
		nbw_str = @sprintf("NBW = %.1f Hz", 𝑓ₛout*1.5/Nfft) #orig. coords: (10, -90); :cl
		push!(plot_𝑓dec,
			cons(:wfrm, WdBw, line=set(style=:solid, color=:red, width=2), label="W"),
			cons(:atext, nbw_str, y=-90, reloffset=set(x=0.05), align=:cl),
		)
	end

	pcoll = push!(cons(:plot_collection, ncolumns=2),
		plot_tdec, plot_𝑓dec, plot_tos, plot_𝑓os,
	)

	return (plot=pcoll, u0=u0, w=w, input=input)
end

end #module DSAudioDemo


#==Auto-run test code:
===============================================================================#
println()
display(@doc(DSAudioDemo)) #Show user how to use DSAudioDemo

function load_demo4_audio_data(m::DSAudioDemo.Modulator)
	srcpath = dirname(pathof(RSDeltaSigmaPort))
	fpath = joinpath(srcpath, "..", "original_source", "delsig", "dsdemo4.mat")
	alldata = DSAudioDemo.MAT.matread(fpath)
	ds = alldata["ds"][:]
	sd = alldata["sd"][:]
	u0 = ds
	u = interp(sd, m.OSR)
	return DSAudioDemo.SrcRaw(u, u0)
end

function playresults(results)
	if isa(results.input, DSAudioDemo.SrcRamp)
		@warn("Will not playback ramp signal.")
		return
	end
	println()
	@info("Listening to ideally sampled/decimated signal...")
		flush(stdout); flush(stderr)
		@show DSAudioDemo.play(results.u0)
	println()
	@info("Listening to modulator output...")
		flush(stdout); flush(stderr)
		@show DSAudioDemo.play(results.w)

	println("\nReplay results with:")
	println("\tplayresults(results)")
	return
end

#Inputs
dsm = DSAudioDemo.Modulator(1, OSR=32, sincorder=2)
#sig = load_demo4_audio_data(dsm)
sig = DSAudioDemo.SrcSine(f=500) #500, 4000, 4200
#sig = DSAudioDemo.SrcRamp()

println("\nUsing modulator:")
@show dsm

println()
@info("Performing ΔΣ audio simulation..."); flush(stdout); flush(stderr)
results = DSAudioDemo.run(dsm, Tsim=3.0, input=sig)
println("\tdone."); flush(stdout); flush(stderr)

println()
@info("Displaying results..."); flush(stdout); flush(stderr)
displaygui(results.plot)

playresults(results)
println()

:END_OF_DEMO
