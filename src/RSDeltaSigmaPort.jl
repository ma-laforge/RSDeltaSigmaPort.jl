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

import Random
import Statistics: mean
import LinearAlgebra
import LinearAlgebra: norm, diagm
import SpecialFunctions: erfinv
import Interpolations
import FFTW: fft, fftshift
import DSP: conv
#import Polynomials #roots, Polynomial
import Printf: @sprintf
using CMDimData
using CMDimData.MDDatasets
using CMDimData.EasyPlot
using CMDimData.EasyPlot: BoundingBox #Use this one to avoid version conflicts with Graphics.
using InspectDR
CMDimData.@includepkg EasyPlotInspect

import ControlSystems
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
include("math.jl")
include("mlfun.jl")
include("mlfun_control.jl")
include("datasets.jl")
include("power.jl")
include("snr.jl")
include("optimize.jl")
include("transferfunctions.jl")
include("statespace.jl")
include("windowing.jl")
include("synthesizeNTF.jl")
include("quantizer.jl")
include("simulateDSM.jl")
include("realizeTF.jl")
include("plotgen.jl")
include("display.jl")

#Data:
export cplxpair, array_round
export _zpk, _zpkdata

#Simple calculations:
export dbm, dbp, dbv, rms
export mapQtoR
export calculateSNR, peakSNR

#Advanced algorithms:
#export fft #From FFTW #Don't export; high risk of collision
export predictSNR

#Transfer functions & windowing:
export evalTF, rmsGain
export synthesizeNTF, realizeNTF
export ds_hann

#State space:
export stuffABCD, scaleABCD, mapABCD

#Simulations:
export simulateSNR, simulateDSM

#Plotting:
export wfrm_stairs, waveform
export plotSpec, plotSpec!
export plotPZ, plotPZ!, plotNTF, plotNTF!
export plotModTransient, plotModSpectrum, plotModSpectrum!
export plotSNR, plotSNR!
export plotStateMaxima
export inlinedisp, saveimage, displaygui

end # module
