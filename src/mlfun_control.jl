#RSDeltaSigmaPort: Implement missing core functionality using ControlSystems.jl
#-------------------------------------------------------------------------------

#=TODO (Goal)
Want to not depend on ControlSystems.jl (slows down compile time significantly).
=#

function _zp2ss(z, p, k)
	H = ControlSystems.zpk(z, p, k)
	return ControlSystems.ss(H)
end

#Last line
