#dsdemo4_audio.jl: A GUI-based demonstration of delta-sigma with sound (Currently GUI-less).

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
using MAT, WAV


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

mutable struct SrcWav <: AbstractSource
	data::Vector #{Float64}
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

function generate(src::SrcWav, m::Modulator; N=10)
	𝑓ₒₛ = m.fsout*m.OSR
#		load('dsdemo4.mat'); #for sd, ds
#		u0 = ds;
#		u = interp(sd,DecFact);
#		N = length(u);
	throw(:INCOMPLETE)
	tidx = (0:N-1); t = tidx/𝑓ₒₛ
	return (t, u, u0)
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
	tplot = cons(:plot, linlin, title="Decimated Signal (dec=$OSR)", legend=true,
		#xyaxes=set(xmin=-120, xmax=0, ymin=0, ymax=120),
		labels=set(xaxis="Time [s]", yaxis="u₀(t)"),
	)

	N = length(u0)
	t = (0:N-1)/𝑓ₛout
	u0w = waveform(t, u0)
	push!(tplot,
		cons(:wfrm, u0w, line=set(style=:solid, color=:blue, width=2), label="u₀"),
	)

	#Plot U(f), from 0 to 𝑓ₛout/2
	𝑓plot = cons(:plot, linlin, title="Spectrum (dec=$OSR)", legend=true,
		xyaxes=set(xmin=0, xmax=𝑓ₛout/2, ymin=-160, ymax=0),
		labels=set(xaxis="Frequency [Hz]", yaxis="U₀(𝑓)"),
	)

	if typeof(input) in [SrcSine, SrcWav]
		N = length(u0)
		local U
		if isa(input, SrcSine)
			U = fft(u0)/(N/4)
		else #SrcWav
			U = fft(u0 .* ds_hann(N))/(N/4)
		end
		𝑓 = range(0, stop=𝑓ₛout, length=N+1); 𝑓 = 𝑓[1:div(N,2)+1]
		UdB = dbv.(U[keys(𝑓)])
		UdBw = waveform(𝑓, UdB)
		push!(𝑓plot,
			cons(:wfrm, UdBw, line=set(style=:solid, color=:blue, width=2), label="U₀"),
		)
	end

	pcoll_in = push!(cons(:plot_collection),
		tplot, 𝑓plot
	)
	displaygui(pcoll_in)

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

	tplot = cons(:plot, linlin, title="Signals @ Modulator Input/Output", legend=true,
		#xyaxes=set(xmin=-120, xmax=0, ymin=0, ymax=120),
		labels=set(xaxis="samples", yaxis="v(t)"),
	)

	_x = collect(0:length(n)-1)
	vw = wfrm_stairs(_x, v[n])
	uw = waveform(_x, u[n])

	push!(tplot,
		cons(:wfrm, vw, line=set(style=:solid, color=:blue, width=2), label="v"),
		cons(:wfrm, uw, line=set(style=:solid, color=:green, width=2), label="u"),
	)

	#Plot V(f), from 0 to Fs/2. Use the middle Nfft points of v
	N = length(v)
	Nfft = min(N, 16*8192)
	n = array_round((N-Nfft)/2+1):array_round((N+Nfft)/2)
	V = fft(v[n] .* ds_hann(Nfft))/(Nfft/4)
	inBin = ceil(Int, Nfft/1000)
	if isa(input, SrcSine)
		inBin = round(Int, input.f/fos*Nfft+1) #Bin of tone
	end
	(𝑓nrm, Vp) = logsmooth(V, inBin)
	#iinf = findall(isinf.(Vp)) #Appears to be infinite values in certain cases
	#finf = 𝑓nrm[iinf]; @show(iinf, finf .* fos)

	𝑓plot = cons(:plot, loglin, title="Spectrum (OSR=$OSR)", legend=true,
		xyaxes=set(xmin=100, xmax=fos/2, ymin=-160, ymax=0),
		labels=set(xaxis="Frequency [Hz]", yaxis="V(𝑓)"),
	)

	Vpw = waveform(𝑓nrm*fos, Vp)
	nbw_str = @sprintf("NBW = %.1f Hz", fos*1.5/Nfft) #orig. coords: (Fs/2, -90); :cr
	push!(𝑓plot,
		cons(:wfrm, Vpw, line=set(style=:solid, color=:blue, width=2), label="V"),
		cons(:atext, nbw_str, y=-90, reloffset=set(x=0.95), align=:cr),
	)

	pcoll_mod = push!(cons(:plot_collection),
		tplot, 𝑓plot
	)
	displaygui(pcoll_mod)

	return
end

end #module ds_audiodemo

display(@doc(DSAudioDemo)) #Show user to use ds_audiodemo

sig = DSAudioDemo.SrcSine(f=500)
#sig = DSAudioDemo.SrcRamp()
#sig = DSAudioDemo.SrcWav()

m = DSAudioDemo.Modulator(1, OSR=32, sincorder=2)
DSAudioDemo.run(m, Tsim=2.0, input=sig)

#Last line
