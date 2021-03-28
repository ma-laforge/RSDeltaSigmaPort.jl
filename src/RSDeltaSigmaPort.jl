#RSDeltaSigmaPort: A port of Richard Schreier's Delta Sigma Toolbox 
#-------------------------------------------------------------------------------

module RSDeltaSigmaPort
#__precompile__(false)
#=Flags:
	-TODO
	-VERIFYME: is this the correct behaviour???
=#

import LinearAlgebra
using Statistics: mean
using LinearAlgebra: norm
using Printf: @sprintf
using CMDimData
using CMDimData.MDDatasets
using CMDimData.EasyPlot
using CMDimData.EasyPlot: BoundingBox #Use this one to avoid version conflicts with Graphics.
using InspectDR
CMDimData.@includepkg EasyPlotInspect

#import ControlSystems
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
include("power.jl")
include("optimize.jl")
include("transferfunctions.jl")
include("synthesizeNTF.jl")
#include("simulateDSM.jl")
include("plotgen.jl")
include("display.jl")

export cplxpair, array_round
export _zpk, _zpkdata
export dbm, dbp, dbv, rms
export evalTF, rmsGain
export synthesizeNTF
export plotPZ, plotPZ!, plotNTF, plotNTF!
export inlinedisp, saveimage, displaygui

end # module
