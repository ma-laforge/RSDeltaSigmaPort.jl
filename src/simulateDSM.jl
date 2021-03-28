#RSDeltaSigma: simulateDSM algorithm
#-------------------------------------------------------------------------------

"""`ds_quantize(y,n)`

Quantize y to:
 - an odd integer in [-n+1, n-1], if n is even, or
 - an even integer in [-n, n], if n is odd.

This definition gives the same step height for both mid-rise
and mid-tread quantizers.
"""
function ds_quantize(y, n::Int)
	if mod(n,2)==0 #mid-rise quantizer
		v = 2*floor(0.5*y)+1
	else #mid-tread quantizer
		v = 2*floor(0.5*(y+1))
	end

	#Limit the output
	for qi in 1:length(n) #Loop for multiple quantizers
		L = n(qi)-1
		i = v(qi,:)>L
		if any(i)
			v(qi,i) = L
		end
		i = v(qi,:)<-L
		if any(i)
			v(qi,i) = -L
		end
	end
	return v
end

#Main algorithm for simulateDSM
function _simulateDSM(u, nq, nlev, x0, order, A, B, C, D1, savestate::Bool, trackmax::Bool)
	if isnan(x0)
		x0 = zeros(order,1)
	end

	N = length(u);
	v = zeros(nq,N);
	y = zeros(nq,N);

	xn = nothing, xmax=nothing
	if savestate #Need to store the state information
		xn = zeros(order,N)
	end
	if trackmax #Need to keep track of the state maxima
		xmax = abs(x0)
	end

	for i=1:N
		y(:,i) = C*x0 + D1*u(:,i);
		v(:,i) = ds_quantize(y(:,i),nlev);
		x0 = A * x0 + B * [u(:,i);v(:,i)];
		if savestate #Save the next state
			xn(:,i) = x0
		end
		if trackmax #Keep track of the state maxima
			xmax = max(abs(x0),xmax)
		end
	end

	return (v,xn,xmax,y)
end


function simulateDSM(u, ntf::ZPKData, nlev=2, x0=NaN, savestate::Bool=false, trackmax::Bool=false)
	nq = length(nlev)
	order = length(ntf.zeros)

	A, B2, C, D2 = zp2ss(ntf.poles,ntf.zeros,-1) #realization of 1/H
	#Transform the realization so that C = [1 0 0 ...]
	Sinv = orth([C' eye(order)])/norm(C); S = inv(Sinv)
	C = C*Sinv
	if C(1)<0
		S = -S
		Sinv = -Sinv
	end
	A = S*A*Sinv; B2 = S*B2; C = [1 zeros(1,order-1)]; #C=C*Sinv
	D2 = 0
	#!!!! Assume stf=1
	B1 = -B2
	D1 = 1
	B = [B1 B2]

	return _simulateDSM(u, nq, nlev, x0, order, A, B, C, D1, savestate, trackmax)
end

function simulateDSM(u, ABCD::Array, nlev=2, x0=NaN)
	nu = size(u,1)
	nq = length(nlev)
	if !(size(ABCD,2) > 2 && size(ABCD,2)==nu+size(ABCD,1))
		msg = "The ABCD argument does not have proper dimensions."
		throw(ArgumentError(msg))
	end
	order = size(ABCD,1)-nq
	A = ABCD(1:order, 1:order);
	B = ABCD(1:order, order+1:order+nu+nq);
	C = ABCD(order+1:order+nq, 1:order);
	D1= ABCD(order+1:order+nq, order+1:order+nu);

	return _simulateDSM(u, nq, nlev, x0, order, A, B, C, D1, savestate, trackmax)
end

"""`(v,xn,xmax,y) = simulateDSM(u,ABCD,nlev=2,x0=0)`
or
`(v,xn,xmax,y) = simulateDSM(u,ntf,nlev=2,x0=0)`

Compute the output of a general delta-sigma modulator with input u,
a structure described by ABCD, an initial state x0 (default zero) and 
a quantizer with a number of levels specified by nlev.
Multiple quantizers are implied by making nlev an array,
and multiple inputs are implied by the number of rows in u.

Alternatively, the modulator may be described by an NTF.
The NTF is zpk object. (The STF is assumed to be 1.)
The structure that is simulated is the block-diagional structure used by
zp2ss.m.
""" simulateDSM

@warn("Implements .m version of simulateDSM. C code might be more efficient.")

#Last line
