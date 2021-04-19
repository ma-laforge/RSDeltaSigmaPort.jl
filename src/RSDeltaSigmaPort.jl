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
   -NEEDSTOOLKIT
	-REQUESTHELP
=#

import Random
import Statistics: mean
import LinearAlgebra
import LinearAlgebra: norm, diagm, eigen
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
include("text.jl")
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
include("calc_spectrum.jl")
include("plot_base.jl")
include("plot_transient.jl")
include("plot_spectrum.jl")
include("plot_NTF.jl")
include("plot_zplane.jl")
include("plot_SNR.jl")
include("plot_state.jl")
include("display.jl")

#Data:
export cplxpair, array_round
export _zpk, _zpkdata
export RealDSM, QuadratureDSM #AbstractDSM?

#Simple calculations:
export dbm, dbp, dbv, rms
export undbm, undbp, undbv
export ds_f1f2
export mapQtoR

#Signal generators:
export default_ftest, default_fband
export genTestTone, genTestTone_quad, genTestTone_sin

#Windowing:
export ds_hann

#Misc. calculations:
export calcSpecInfo
export calcSNRInfo, calculateSNR, peakSNR, predictSNR

#Transfer functions:
export evalTF, rmsGain
export synthesizeNTF, realizeNTF

#State space:
export stuffABCD, scaleABCD, mapABCD

#Simulations:
export simulateDSM, simulateSNR

#Plotting:
export wfrm_stairs, waveform
export plotSpec, plotSpec!
export plotPZ, plotPZ!, plotNTF, plotNTF!
export plotModTransient, plotModSpectrum, plotModSpectrum!
export plotExampleSpectrum
export plotSNR, plotSNR!
export documentNTF
export plotStateMaxima

#Displaying plots:
export inlinedisp, saveimage, displaygui

#Display/text:
export ds_orderString, str_modulatortype

end # module
