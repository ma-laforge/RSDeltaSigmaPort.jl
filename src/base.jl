#RSDeltaSigma: Base functionnality and types
#-------------------------------------------------------------------------------

#=Info
  - Circumvent using ControlSystems.jl for now by specifying custom zpk()
    functions & data structures.
=#


#==Type definitions
===============================================================================#
const IndexRange = UnitRange{Int64}

const MISSINGARG = ArgumentError("Missing Argument")


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

abstract type AbstractDSM; end

"""`RealDSM(;order=3, OSR=16, M=1, f0=0, form=:CRFB, Hinf=1.5, opt=0)`

Describe a real-valued (not quadrature) ΔΣ modulator.

# Inputs
 - `order`: order of the modulator
 - `OSR`: oversampling ratio
 - `M`: Number of quantizer steps (`M=nlev-1`)
 - `f0`: center frequency (1->fs)
 - `form`:
 - `H_inf`: maximum NTF gain
 - `opt`: flag for optimized zeros\r
 - `--> opt=0` -> not optimized,
 - `--> opt=1` -> optimized,
 - `--> opt=2` -> optimized with at least one zero at band-center,
 - `--> opt=3` -> optimized zeros (Requires MATLAB6 and Optimization Toolbox),
 - `--> opt=[z]` -> zero locations in complex form (TODO?)
"""
struct RealDSM <: AbstractDSM
	order::Int
	OSR::Int
	M::Int #Num steps
	f0::Float64
	form::Symbol
	Hinf::Float64
	opt::Int
end
RealDSM(;order=3, OSR=16, M=1, f0=0, form=:CRFB, Hinf=1.5, opt=0) =
	RealDSM(order, OSR, M, f0, form, Hinf, opt)

struct QuadratureDSM <: AbstractDSM
	order::Int
	OSR::Int
	M::Int #Num steps
	f0::Float64
	form::Symbol
	NG::Int
	ING::Int
end
QuadratureDSM(;order=4, OSR=32, M=1, f0=1/16, form=:PFB, NG=-50, ING=-10) =
	QuadratureDSM(order, OSR, M, f0, form, NG, ING)


#==Accessors
===============================================================================#
isquadrature(::RealDSM) = false
isquadrature(::QuadratureDSM) = true


#==Helper functions
===============================================================================#
function orderinfo(order::Int)
	order2 = div(order, 2)
	Δodd = mod(order, 2)
	isodd = order > 2*order2
	return (order2, isodd, Δodd)
end
orderinfo(dsm::AbstractDSM) = orderinfo(dsm.order)


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


#==Basic frequency calculations/defaults
===============================================================================#
"""`(f1, f2) = ds_f1f2(;OSR=64, f0=0, iscomplex=0)`

Deprecated. Use `default_fband()` instead.

# See also:
[`default_fband`](@ref), [`default_ftest`](@ref)
"""
function ds_f1f2(;OSR::Int=64, f0::Float64=0, iscomplex::Bool=false)
	local f1, f2
	if iscomplex
		f1 = f0-0.5/OSR
		f2 = f0+0.5/OSR
	else
		if f0>0.25/OSR
			f1 = f0-0.25/OSR
			f2 = f0+0.25/OSR
		else
			f1 = 0
			f2 = 0.5/OSR
		end
	end
	return (f1, f2)
end

"""`[f1, f2] = default_fband(dsm)`

Compute a default application band (normalized frequency limits) for the given
ΔΣ modulator.

# See also:
[`default_ftest`](@ref)
""" default_fband

default_fband(OSR::Int; f0::Float64=0.0, quadrature::Bool=false) =
	collect(ds_f1f2(OSR=OSR, f0=f0, iscomplex=quadrature)) #Vector form
default_fband(dsm::RealDSM) = default_fband(dsm.OSR, f0=dsm.f0, quadrature=false)
default_fband(dsm::QuadratureDSM) = default_fband(dsm.OSR, f0=dsm.f0, quadrature=false)

"""`default_ftest(dsm)`

Compute a default test (normalized) frequency for the given ΔΣ modulator.

# See also:
[`default_fband`](@ref)
""" default_ftest

function default_ftest(OSR::Int; f0::Float64=0, quadrature::Bool=false)
	if quadrature
		return ftest = 0.15/OSR
	else
		return (0==f0) ? 0.15/OSR : f0 + 0.08/OSR
	end
end
default_ftest(dsm::RealDSM) = default_ftest(dsm.OSR, f0=dsm.f0, quadrature=false)
default_ftest(dsm::QuadratureDSM) = default_ftest(dsm.OSR, f0=dsm.f0, quadrature=true)

#Last line
