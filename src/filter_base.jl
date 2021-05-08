#RSDeltaSigma: Base functions to assist with filtering operations
#-------------------------------------------------------------------------------

#==Types
===============================================================================#
"""`CoeffCSD`

A structure containing filter coefficents and their canonical-signed digit (CSD)
representation.
"""
struct CoeffCSD
	val
	csd
end


#==Accessors
===============================================================================#
csdsize2(f::Vector{CoeffCSD}) = [size(fi.csd,2) for fi in f] #VERIFYME

Base.values(f::Vector{CoeffCSD}) = [fi.val for fi in f]


#==Functions
===============================================================================#
"""`bquantize(x, nsd=3, abstol=eps, reltol=10*eps)`

Bidirectionally quantize a n by 1 vector x to nsd signed digits.

Terminate early if the error is less than the specified tolerances.

# Output
y is a CoeffCSD[] with the following fields:
 - `y[i].val` is the quantized value in floating-point form
 - `y[i].csd` is a 2-by-nsd (or less) matrix containing
   the powers of two (first row) and their signs (second row).

# See also:
[`bunquantize`](@ref)
"""
function bquantize(x, nsd::Int=3, abstol::Float64=eps(1.0), reltol::Float64=10*eps(1.0))
	n = length(x)
	q = zeros(2*n,nsd)
	y = Array{CoeffCSD}(undef, n)
	offset = -log2(0.75)

	for i in 1:n
		xp = x[i]
		val = 0
		csd = zeros(2,nsd)
		for j in 1:nsd
			error = abs(val-x[i])
			if error <= abstol || error <= abs(x[i])*reltol
				break
			end
			p = floor( log2(abs(xp)) + offset )
			p2 = 2^p
			sx = sign(xp)
			xp = xp - sx*p2
			val = val + sx*p2
			csd[1:2, j] = [p; sx]
		end
		y[i] = CoeffCSD(val, csd)
	end
	return y
end

"""`bunquantize(q)`

Calculate the value corresponding to a bidirectionally quantized quantity.
q is a 2n by m matrix containing the powers of 2 and their signs for each
quantized value.

# See also:
[`bquantize`](@ref)
"""
function bunquantize(q)
	q = 1.0*q #Ensure we have floating point values
	n = size(q,1)/2
	y = zeros(array_round(n),1)
	signs = 2:2:array_round(2*n)
	powers = signs .- 1
	for i in 1:size(q,2)
		y = y + 2 .^ q[powers,i] .* q[signs,i]
	end
	return y
end

function convertForm(f,q)
	n = length(f)
	F = CoeffCSD[]
	for i in 1:n
		val = f[i]
		csdrep = q[2*i-1:2*i,:]
		csd = csdrep[:,findall(csdrep[1,:] .!=0)]
		push!(F, CoeffCSD(val, csd))
	end
	return F
end

#Last line
