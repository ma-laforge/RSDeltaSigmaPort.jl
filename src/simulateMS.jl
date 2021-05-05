#RSDeltaSigma: simulateMS() algorithm
#-------------------------------------------------------------------------------


#==simulateMS()
===============================================================================#
"""`simulateMS(v, M=16, mtf=zpk(1,0,1,1), d=0, dw=[1..], sx0=[0..])`

Simulate the vector-based Mismatch-Shaping for a multi-element DAC:
 - simresult = (sv, sx, sigma_se, max_sx, max_sy)

# Inputs:
 - `v`: A vector of the digital input values. v is in (-M:2:M) if dw=[1..]
   otherwise v is in [-sum(dw),sum(dw)].
 - `M`: The number of elements.
 - `mtf`: The mismatch-shaping transfer function, given in zero-pole form.
 - `d`: Dither uniformly distributed in [-d,d] is added to sy.
 - `dw`: A vector of dac element weights
 - `sx0`: A matrix whose columns are the initial states of the ESL.

# Outputs:
 - `sv`: An MxN matrix whose columns are the selection vectors.
 - `sx`: An orderxM matrix containing the final state of the ESL.
 - `sigma_se`: The rms value of the selection error.
 - `max_sx`: The maximum absolute value of the state for all modulators.
 - `max_sy`: The maximum absolute value of the input to the VQ.

# Notes:
For a brief description of the theory of mismatch-shaping DACs, see
R. Schreier and B. Zhang "Noise-shaped multibit D/A convertor employing
unit elements", Electronics Letters, vol. 31, no. 20, pp. 1712-1713,
Sept. 28 1995.
"""
function simulateMS(v; M::Int=16, mtf=[], d::Float64=0.0, dw::Vector=[], sx0=[])
	local mtf_z; local mtf_p; local mtf_k

	if isa(mtf, Array) && isempty(mtf)
		mtf_z = [1]; mtf_p = [0]
	elseif isa(mtf, ZPKData)
		(mtf_z, mtf_p, mtf_k) = _zpkdata(mtf)
	else
		_type = typeof(mtf)
		throw("mtf type not supported: $_type.")
	end

	order = length(mtf_p)
	if isempty(sx0)
		sx0 = zeros(order,M)
	end
	if isempty(dw)
		dw = ones(M)
	end

	#B/A = MTF-1
	num = poly(mtf_z)
	den = poly(mtf_p)
	A = real.(-den[2:order+1])
	B = real.(num[2:order+1]+A)
	A = A'; B = B' #Need row vectors
	N = length(v)
	sv = zeros(M,N)

	sx = sx0
	max_sx = maximum(maximum(abs.(sx)))
	max_sy = 0
	sum_se2 = 0

	for i in 1:N
		#Compute the sy vector.
		sy = B*sx
		#Normalize sy for a minimum value of zero.
		sy = sy .- minimum(sy)
		syd = (sy + d*(2*rand(1,M) .- 1))[:]
		#Pick the elements that have the largest desired_usage (sy) values.
		sv[:,i] = selectElement(v[i], syd, dw)
		#Compute the selection error.
		se = sv[:,[i]]' - sy
		#Compute the new sx matrix
		sxn = A*sx + se
		sx = [ sxn; sx[1:order-1,:]]
		#Keep track of some statistics.
		sum_se2 = sum_se2 + sum((se .- mean(se)) .^ 2)
		max_sx = maximum([max_sx abs.(sxn)])
		max_sy = maximum([max_sy abs.(sy)])
	end
	sigma_se = sqrt(sum_se2/(M*N))

	simresult = (sv=sv, sx=sx, sigma_se=sigma_se, max_sx=max_sx, max_sy=max_sy)
	return simresult
end

@warn("Implements .m version of simulateMS(). C code might be more efficient/validated.")

#Last line
