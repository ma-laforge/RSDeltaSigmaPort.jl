#RSDeltaSigmaPort: utilities for transfer functions
#-------------------------------------------------------------------------------

"""`y = evalRPoly(roots,x,k=1)`
Compute the value of a polynomial which is given in terms of its roots.
"""
function evalRPoly(roots,x,k=1)
	y = k * ones(size(x))
	finiteroots = findall(isfinite.(roots))
	roots = roots[finiteroots]
	for i in 1:length(roots)
		y = y.*(x .- roots[i])
	end
	return y
end

"""`h = evalTF(tf,z)`
Evaluates the rational function described by the struct tf
at the point(s) given in the z vector.
TF must be either a zpk object or a struct containing
  form		'zp' or 'coeff'
  zeros,poles,k	if form=='zp'
  num,den		if form=='coeff'

In Matlab 5, the ss/freqresp() function does nearly the same thing.
"""
function evalTF(tf::ZPKData,z)
	h = tf.k * evalRPoly(tf.z,z) ./ evalRPoly(tf.p,z)
	return h
end

function evalTF(tf,z)
	throw(ErrorException("Transfer function of this type not currently supported."))
#elseif any(strcmp(fieldnames(tf),'form'))
#    if strcmp(tf.form,'zp')
#        h = tf.k * evalRPoly(tf.zeros,z) ./ evalRPoly(tf.poles,z);
#    elseif strcmp(tf.form,'coeff')
#        h = polyval(tf.num,z) ./ polyval(tf.den,z);
#    else
#        fprintf(1,'%s: Unknown form: %s\n', mfilename, tf.form);
#    end
#else	% Assume zp form
#    h = tf.k * evalRPoly(tf.zeros,z) ./ evalRPoly(tf.poles,z);
#end
end


#==Calculators
===============================================================================#

"""`g = rmsGain(H,f1,f2,N=100)`

Compute the root mean-square gain of the discrete-time
tf H in the frequency band (f1,f2)
"""
function rmsGain(H,f1,f2,N::Int=100)
	w = collect(range(2*pi*f1, stop=2*pi*f2, length=N))
	g = norm( evalTF(H, exp.(j*w)) ) / sqrt(N)
	return g
end


#Last line
