#RSDeltaSigma: simulateHBF() algorithm
#-------------------------------------------------------------------------------

"""`y = simulateHBF(x, f1, f2, mode=0)`

Simulate a Saramaki half-band filter.

The `f1` and `f2` vectors contain coefficients for the structure.
(`f1` and `f2` can also be struct arrays like those returned from designHBF.m
i.e. struct arrays whose .val fields contain the coeffiecients.)

The `mode` flag determines whether the input is filtered, interpolated,
or decimated according to the following table:
 - `mode = 0`: Plain filtering, no interpolation or decimation.
 - `mode = 1`: The input is interpolated.
 - `mode = 2`: The output is decimated, even samples are taken.
 - `mode = 3`: The output is decimated, odd samples are taken.
"""
function simulateHBF(x, f1, f2, mode::Int=0)
	if isa(f1, Number) && isnan(f1)
		(f1, f2) = exampleHBF(4)
	end
	if isa(f1, Vector{CoeffCSD})
		f1 = values(f1)
	end
	if isa(f2, Vector{CoeffCSD})
		f2 = values(f2)
	end
	x = x[:]
	f1 = f1[:]; f2 = f2[:] #Redundant unless `val`s are not always simple numbers
	n1 = length(f1); n2 = length(f2)
	f2imp = vcat(f2[end:-1:1], f2)
	if 0 == mode #Insert zeros
		f2imp = vcat(f2imp, zeros(2*n2))
		f2imp = f2imp[1:4*n2-1]
	end
	F2 = _tf(f2imp, vcat(1, zeros(length(f2imp)-1)), 1)

	if 0 == mode #Plain
		(up,) = _lsim(F2,x)
		y = 0.5*delay(x,2*n2-1) + f1[1]*up
		for i in 2:n1
			(up,) = _lsim(F2,up)
			(up,) = _lsim(F2,up)
			y = f1[i]*up + delay(y, 4*n2-2)
		end
	elseif 1 == mode #Interpolating
		up = zeros(size(x))
		nz = 2*n2-1
		for i in n1:-1:1
			if i==1
				(up,) = _lsim(F2,up+f1[i]*x)
				nz = n2-1
			else
				(up,) = _lsim(F2,up+f1[i]*x)
				(up,) = _lsim(F2,up)
			end
			x = delay(x,nz)
		end
		y = [2*up'; x']; y = y[:] #Interleave the upper and lower streams

	elseif mode in [2, 3] #Decimating
		x = x[1:2*div(end,2)]
		if 3 == mode
			y = 0.5*x[1:2:end]
			up = x[2:2:end]
			nz = n2-1
		else
			y = 0.5*x[2:2:end]
			up = x[1:2:end]
			nz = n2
		end
		for i in 1:n1
			if i==1
				(up,) = _lsim(F2,up)
			else
				(up,) = _lsim(F2,up)
				(up,) = _lsim(F2,up)
			end
			y = f1[i]*up + delay(y,nz)
			nz = 2*n2-1
		end
	else
		throw("$mode is not a valid value for the mode variable.")
	end
	return y
end

#Last line
