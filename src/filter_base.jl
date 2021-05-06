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
