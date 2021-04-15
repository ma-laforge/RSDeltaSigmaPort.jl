#RSDeltaSigma: Base functionnality and types
#-------------------------------------------------------------------------------

#=Info
  - Circumvent using ControlSystems.jl for now by specifying custom zpk()
    functions & data structures.
=#


#==Type definitions
===============================================================================#
const IndexRange = UnitRange{Int64}

#Type used to dispatch on a symbol & minimize namespace pollution:
struct DS{Symbol}; end; #Dispatchable symbol
DS(v::Symbol) = DS{v}()

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


#==Helper functions
===============================================================================#
function orderinfo(order::Int)
	order2 = div(order, 2)
	Δodd = mod(order, 2)
	isodd = order > 2*order2
	return (order2, isodd, Δodd)
end


#Sort roots such that real values show up first
function realfirst(a)
	epstol = 100

	#Index of real roots:
	Δmin = epstol*eps.(abs.(a)) #Threshold for numerical error
	ireal = abs.(imag.(a)) .< Δmin

	areal = real.(a[ireal]) #Clean up real roots
	acplx = a[.!ireal] #complex pairs
	return vcat(sort(areal), acplx) #Sort roots
end


#==Converters to make input data uniform
===============================================================================#
"""`conv2seriesmatrix2D(x)`

Ensure we have a 2D matrix representing series data (Vector->2×N Array)
"""
conv2seriesmatrix2D(x::T) where T =
	throw(ErrorException("Cannot convert $T to a \"series data matrix\""))
conv2seriesmatrix2D(x::AbstractVector) = collect(x')
conv2seriesmatrix2D(x::Vector) = collect(x') #Otherwise Array traps it.
conv2seriesmatrix2D(x::Array) = x #Assume format is ok.

#Last line
