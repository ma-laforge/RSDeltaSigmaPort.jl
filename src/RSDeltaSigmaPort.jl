#RSDeltaSigmaPort: A port of Richard Schreier's Delta Sigma Toolbox 
#-------------------------------------------------------------------------------
#__precompile__(false) #Bit faster to develop when using ControlSystems.jl
#=Conventions
 - Function name `_stairs__()`
   - Leading `_`: Private function (not meant to be called externally from this
     module).
   - Trailing `__`: Implementation for a specific function (do not call unless
     from within that parent function).
=#

"""`RSDeltaSigmaPort`: A port of Richard Schreier's Delta Sigma Toolbox

# Sample usage
```julia-repl
julia> using RSDeltaSigmaPort #Will take a while to load, compile, etc...
julia> import RSDeltaSigmaPort: @runsample

julia> @runsample("dsdemo1.jl")
julia> @runsample("dsdemo2.jl")
julia> @runsample("dsdemo3.jl")
julia> @runsample("dsdemo4_audio.jl")
julia> @runsample("dsdemo5.jl")
julia> @runsample("dsdemo6.jl")
julia> @runsample("dsexample1.jl")
julia> @runsample("dsexample2.jl")
julia> @runsample("demoLPandBP.jl")
```

# See also:
[`simulateDSM`](@ref), [`simulateMS`](@ref), [`simulateSNR`](@ref), [`simulateHBF`](@ref) | 
[`synthesizeNTF`](@ref), [`realizeNTF`](@ref), [`realizeNTF_ct`](@ref) | 
[`calculateSNR`](@ref), [`peakSNR`](@ref), [`predictSNR`](@ref) | 
[`calculateTF`](@ref), [`evalTF`](@ref), [`evalTFP`](@ref) | 
[`stuffABCD`](@ref), [`scaleABCD`](@ref), [`mapABCD`](@ref), [`partitionABCD`](@ref) | 
[`mapCtoD`](@ref), [`mapQtoR`](@ref) | 
[`exampleHBF`](@ref), [`designHBF`](@ref), [`simulateHBF`](@ref) | 
[`pulse`](@ref), [`impL1`](@ref) | 
[`lollipop`](@ref), [`logsmooth`](@ref) | 
[`documentNTF`](@ref), [`plotExampleSpectrum`](@ref)
"""
module RSDeltaSigmaPort

const rootpath = realpath(joinpath(@__DIR__, ".."))

#=Flags:
	-TODO
	-VERIFYME: is this the correct behaviour???
   -VERSIONSWITCH: Algorithm originally depends on software version (Needs to be checked)
   -NEEDSTOOLKIT
	-REQUESTHELP
=#

import Random
import Statistics: mean
import LinearAlgebra
import LinearAlgebra: norm, diagm, eigen, cond
import SpecialFunctions: erfinv
import Interpolations
import FFTW: fft, fftshift
import DSP: conv, filt, remez
import Optim
import Optim: optimize, GoldenSection, IPNewton, GradientDescent, LBFGS
import Optim: Fminbox, TwiceDifferentiableConstraints
import Polynomials
import Polynomials: Polynomial
import Printf: @sprintf
using CMDimData
using CMDimData.MDDatasets
using CMDimData.EasyPlot
using CMDimData.EasyPlot.Colors #Should be ok with whatever version it needs
import CMDimData.EasyPlot: BoundingBox #Use this one to avoid version conflicts with Graphics.
using InspectDR
CMDimData.@includepkg EasyPlotInspect

import ControlSystems
import ControlSystems: isdiscrete
#using ControlSystems: zpk, zpkdata

function throw_unreachable()
	msg = "Code expected to be unreachable with exception handling"
	error(msg)
end

function throw_notimplemented()
	msg = "Feature not implemented"
	error(msg)
end

const j = im
include("base.jl")
include("arrays.jl")
include("mlfun.jl")
include("mlfun_control.jl")
include("text.jl")
include("datasets.jl")
include("signals_time.jl")
include("power.jl")
include("timedomain.jl")
include("windowing.jl")
include("snr.jl")
include("optimize.jl")
include("transferfunctions.jl")
include("statespace.jl")
include("mapCtoD.jl")
include("synthesizeNTF.jl")
include("quantizer.jl")
include("simulate_base.jl")
include("simulateDSM.jl")
include("simulateMS.jl")
include("realizeNTF.jl")
include("realizeNTF_ct.jl")
include("filter_base.jl")
include("exampleHBF.jl")
include("designHBF.jl")
include("designLCBP.jl")
include("simulateHBF.jl")
include("calc_spectrum.jl")
include("plot_base.jl")
include("plot_transient.jl")
include("plot_spectrum.jl")
include("plot_NTF.jl")
include("plot_zplane.jl")
include("plot_SNR.jl")
include("plot_state.jl")
include("display.jl")


#==Convenience macros
===============================================================================#
macro runsample(filename::String)
	fnstr = "$filename"
	fpath = realpath(joinpath(rootpath, "sample", fnstr))
	@info("Running RSDeltaSigmaPort sample:\n$fpath...")
	m = quote
		include($fpath)
	end
	return esc(m) #esc: Evaluate in calling module
end


#==Exported interface
===============================================================================#

#Simple calculations:
export dbm, dbp, dbv, rms
export undbm, undbp, undbv

#Compatibility:
#WARN: Might eventually collide with names from other packages:
#NOTE: Some use "_" prefix to avoid collisions with similar functions in other packages
export cplxpair
export _ss, _zpk, _zpkdata
export _zp2ss, _zp2tf, _freqz, _tf, _minreal, _impulse
export interp1_lin, interp1_cubic #Specialized names
export interp
export squeeze
export eye, orth, eig

#Data/arrays:
export AbstractDSM, RealDSM, QuadratureDSM
export padl, padr, padt, padb
export array_round #Converts a floating-point value to an (Int) index usable by Julia arrays.
export mapQtoR

#Accessors
export isquadrature

#Signal generators:
export ds_f1f2, default_ftest, default_fband
export genTestTone, genTestTone_quad, genTestTone_sin
export pulse, impL1

#Windowing/time domain:
export ds_hann, applywnd
export ds_therm
export sinc_decimate

#Misc. calculations:
export calcSpecInfo
export calcSNRInfo, calculateSNR, peakSNR, predictSNR

#Transfer functions:
export synthesizeNTF, realizeNTF, realizeNTF_ct
export calculateTF, evalTF, evalTFP
export rmsGain, cancelPZ

#State space:
export stuffABCD, scaleABCD, mapABCD, partitionABCD
export mapCtoD

#Filters:
export exampleHBF, designHBF, designLCBP

#Simulations:
export simulateDSM, simulateMS
export simulateSNR, simulateHBF

#Plotting:
export waveform, wfrm_stairs, lollipop, logsmooth
export plotPZ, plotPZ!, plotNTF, plotNTF!
export documentNTF
export plotSNR, plotSNR!, plotStateMaxima, plotUsage
export plotModTransient, plotLollipop
export plotSpec, plotSpec!, plotModSpectrum, plotModSpectrum!
export plotSpectrum, plotSpectrum!, plotExampleSpectrum

#Displaying plots:
export inlinedisp, saveimage, displaygui

#Display/text:
export ds_orderString, str_modulatortype

end # module
