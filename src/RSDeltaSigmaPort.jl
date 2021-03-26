module RSDeltaSigmaPort
#__precompile__(false)
#=Flags:
	-TODO
	-VERIFYME: is this the correct behaviour???
=#

using Graphics: BoundingBox
using Statistics: mean
using LinearAlgebra: norm
using Printf: @sprintf
using CMDimData
using CMDimData.MDDatasets
using CMDimData.EasyPlot

#using InspectDR
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
include("helperfunc.jl")
include("transferfunctions.jl")
include("synthesizeNTF.jl")
include("plotgen.jl")

export dbm, dbp, dbv, rms
export evalTF, rmsGain
export synthesizeNTF
export plotPZ, plotPZ!, plotNTF, plotNTF!

end # module
