#RSDeltaSigmaPort: A port of Richard Schreier's Delta Sigma Toolbox 
#-------------------------------------------------------------------------------
__precompile__(false) #Bit faster to develop when using ControlSystems.jl
#=Conventions
 - Function name `_stairs__()`
   - Leading `_`: Private function (not meant to be called externally from this
     module).
   - Trailing `__`: Implementation for a specific function (do not call unless
     from within that parent function).
=#


module RSDeltaSigmaPort
#=Flags:
	-TODO
	-VERIFYME: is this the correct behaviour???
=#

import LinearAlgebra
import FFTW
using Statistics: mean
using LinearAlgebra: norm
using FFTW: fft
using Printf: @sprintf
using CMDimData
using CMDimData.MDDatasets
using CMDimData.EasyPlot
using CMDimData.EasyPlot: BoundingBox #Use this one to avoid version conflicts with Graphics.
using InspectDR
CMDimData.@includepkg EasyPlotInspect

import ControlSystems
import ControlSystems.Polynomials #roots, Polynomial
#using ControlSystems: zpk, zpkdata

"Availability of `fmincon` function from the Optimization Toolbox"
const FMINCON_AVAIL = false

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
include("mlfun.jl")
include("mlfun_control.jl")
include("datasets.jl")
include("power.jl")
include("snr.jl")
include("optimize.jl")
include("transferfunctions.jl")
include("windowing.jl")
include("synthesizeNTF.jl")
include("quantizer.jl")
include("simulateDSM.jl")
include("plotgen.jl")
include("display.jl")

#Data:
export cplxpair, array_round
export _zpk, _zpkdata

#Simple calculations:
export dbm, dbp, dbv, rms
export calculateSNR

#Advanced algorithms:
export fft #From FFTW

#Transfer functions & windowing:
export evalTF, rmsGain
export synthesizeNTF
export ds_hann

#Simulations:
export simulateDSM

#Plotting:
export wfrm_stairs, waveform
export plotSpec, plotSpec!
export plotPZ, plotPZ!, plotNTF, plotNTF!
export plotModTransient, plotModSpectrum
export inlinedisp, saveimage, displaygui

end # module
