#RSDeltaSigma: Base functionnality and types
#-------------------------------------------------------------------------------

#=Info
  - Circumvent using ControlSystems.jl for now by specifying custom zpk()
    functions & data structures.
=#


#==Type definitions
===============================================================================#
const IndexRange = UnitRange{Int64}

"""`ZPKData`

Custom structure to store zpk data.
"""
mutable struct ZPKData
	z::Union{Float64, Vector{Float64}, Vector{Complex{Float64}}}
	p::Union{Float64, Vector{Float64}, Vector{Complex{Float64}}}
	k::Real
	Ts::Real
end

function _zpk(z, p, k::Number, Ts::Number)
	ZPKData(z, p, k, Ts)
end

_zpkdata(d::ZPKData) = (d.z, d.p, d.k)


#==Converters to make input data uniform
===============================================================================#
"""`conv2seriesmatrix2D(x)`

Ensure we have a 2D matrix representing series data (Vector->2×N Array)
"""
conv2seriesmatrix2D(x::T) where T =
	throw(ErrorException("Cannot convert $T to a \"series data matrix\""))
conv2seriesmatrix2D(x::Vector) = collect(x')
conv2seriesmatrix2D(x::Array) = x #Assume format is ok.

#Last line
