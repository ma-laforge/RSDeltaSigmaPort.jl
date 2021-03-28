#RSDeltaSigmaPort: Implement missing core functionality
#-------------------------------------------------------------------------------

"""Emulates rounding of array indices"""
array_round(i) = round(Int, i) #VERIFYME

function cplxpair(a)
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

#Last line
