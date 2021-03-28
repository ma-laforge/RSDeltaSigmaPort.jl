#RSDeltaSigma: Base functionnality and types
#-------------------------------------------------------------------------------

#=Info
  - Circumvent using ControlSystems.jl for now by specifying custom zpk()
    functions & data structures.
=#

#Custom structure to store zpk data.
mutable struct ZPKData
	z::Union{Float64, Vector{Float64}, Vector{Complex{Float64}}}
	p::Union{Float64, Vector{Float64}, Vector{Complex{Float64}}}
	k::Real
	Ts::Real
end

function _zpk(z, p, k::Number, Ts::Number)
	ZPKData(z, p, k, Ts)
end

_zpkdata(d::ZPKData) = (d.z, d.p, dpk)


#Last line
