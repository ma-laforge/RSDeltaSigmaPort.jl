#RSDeltaSigmaPort: Implement missing core functionality
#-------------------------------------------------------------------------------

"""Emulates rounding of array indices"""
array_round(i) = round(Int, i) #VERIFYME

function rat(v::Real, tol::Float64)
	r=rationalize(Float64(v), tol=tol)
	return (r.num, r.den)
end

function squeeze(a)
	sz = collect(size(a))
	for i in reverse(keys(sz))
		if sz[i]<2
			a = dropdims(a, dims=i)
		end
	end
	return a
end

#Sort complex poles
function _cplxpair(a)
	epstol = 100

	#Sorts real first, then imaginary:
	function isrealless(a,b)
		Δmin = epstol*eps(abs(a))
		if abs(real(a) - real(b)) < Δmin
			return imag(a) < imag(b)
		else
			return real(a) < real(b)
		end
	end

	a = sort(a, lt=isrealless)
	i = 1
	while i < length(a)
		v = a[i]
		Δmin = epstol*eps(abs(v))
		if abs(imag(v)) < Δmin #Really a real number
			a[i] = real(v)
			i+=1
			break
		elseif i < length(a)
			if abs(v-conj(a[i+1])) < Δmin
				i+=2 #Skip both cplx conj pairs.
				break
			else
				throw(ErrorException("cplxpair failed."))
			end
		else
			throw(ErrorException("cplxpair failed."))
		end
		
	end
	return a
end

"""`cplxpair(a)`

Sort roots in array a such that complex pairs are in ascending order of their
real parts, then by their imaginary part.  In the output, the real roots
re-positionned after the complex roots.
"""
function cplxpair(a)
	epstol = 100

	#Index of real roots:
	Δmin = epstol*eps.(abs.(a)) #Threshold for numerical error
	ireal = abs.(imag.(a)) .< Δmin

	areal = real.(a[ireal]) #Clean up real roots
	acplx = a[.!ireal] #complex pairs
	return vcat(_cplxpair(acplx), sort(areal)) #Sort roots
end

#Custom function to ensure complex conjugate pairs match exactly.
function cleancplxpair(l)
	epstol = 100

	if length(l) < 2
		return l
	end

	l = cplxpair(l) #Sort to group complex conjugates together
	result = similar(l)
	if length(l)>0
		#Copy last value in case not part of conjugate pair
		#(ignored by following algorithm)
		result[end] = l[end]
	end

	i = 1
	while i < length(l)
		a, b = l[i], l[i+1]
		ΔRmin = epstol*eps(real(a))
		ΔXmin = epstol*eps(imag(a))

		result[i] = a
		if abs(real(a) - real(b)) < ΔRmin && abs(imag(a) + imag(b)) < ΔXmin
			a = (a+conj(b))/2 #Take mean & re-write
			result[i] = a
			result[i+1] = conj(a)
			i+=1
		end
		i+=1
	end

	return result
end

#TODO: Will algorithms work without collect (possibly faster?):
eye(n::Int) = collect(LinearAlgebra.I(n))

#It sounds like this should be implementation of orth:
#orth(A) = LinearAlgebra.qr(A).Q

#But scipy seems to do it this way:
orth(A) = LinearAlgebra.svd(A).U

#Convert array of roots into polynomial:
function poly(rA)
	result = [1]
	for r in rA
		result = conv(result, [1, -r])
	end
	return result
end

eig(a) = eigen(a).values

#Implement interp1 using linear interpolation:
function interp1_lin(xA, yA, xv)
	ILib = Interpolations
	itp = ILib.interpolate((x,), y, ILib.Gridded(ILib.Linear()))
		return itp(xv)
end

#Implement interp1 using cubic interpolation:
function interp1_cubic(xA, yA, xv)
	ILib = Interpolations
   itp = ILib.interpolate(yA, ILib.BSpline(ILib.Cubic(ILib.Natural())), ILib.OnGrid())
		intf = ILib.scale(itp, xA)
	return intf(xv)
end

#FIXME/VERIFYME
#Algoritm assembled together by MALaforge, hoping it will be sufficient.
#Probably a good idea to have an expert implement proper function.
function interp(v, OSR::Int)
	thresh = .01
	Δs = OSR*ceil(Int, 1/(thresh*π))
	x = -Δs:1:Δs
	v_0 = zeros(length(v)*OSR)
	for i in 0:length(v)-1
		v_0[1+OSR*i] = v[1+i]
	end
	_filt = sinc.(x/OSR)
	result = conv(_filt, v_0)
	result = result[Δs .+ (1:length(v_0))]
	return result
end

#Evaluate polynomial p at value x
function polyval(p, x)
	#NOTE: [1, 2, 3] => 1 + 2x¹ + 3x²
	polyfn = Polynomial(p)
	return polyfn.(x) #Evaluates polynomial
end

#_roots so as to not confuse it with Polynomials.roots
#NOTE: Polynomial expects the coefficients in opposite order:
_roots(p) = Polynomials.roots(Polynomial(reverse(p)))

#Last line
