module RSDeltaSigmaPort
#__precompile__(false)

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

include("base.jl")
include("mlfun.jl")
include("conversion.jl")
include("helperfunc.jl")
include("transferfunctions.jl")
include("synthesizeNTF.jl")

export dbm, dbp, dbv
export evalTF
export synthesizeNTF

end # module
