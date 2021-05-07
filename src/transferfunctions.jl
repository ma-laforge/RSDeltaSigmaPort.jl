#RSDeltaSigmaPort: utilities for transfer functions
#-------------------------------------------------------------------------------

"""`y = evalRPoly(roots,x,k=1)`
Compute the value of a polynomial which is given in terms of its roots.
"""
function evalRPoly(roots,x::AbstractArray,k=1)
	y = k * ones(size(x))
	finiteroots = findall(isfinite.(roots))
	roots = roots[finiteroots]
	for i in 1:length(roots)
		y = y.*(x .- roots[i])
	end
	return y
end
evalRPoly(roots,x,k=1) = evalRPoly(roots, [x], k)[1]

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


"""`H=evalTFP(Hs, Hz, f)`

Compute the value of a transfer function product Hs*Hz at a frequency f,
where Hs is a cts-time TF and Hz is a discrete-time TF.
Both Hs and Hz are SISO zpk objects.
This function attempts to cancel poles in Hs with zeros in Hz.
"""
function evalTFP(Hs,Hz,f::AbstractArray)
	(szeros, spoles, sk) = _zpkdata(Hs, 1)
	(zzeros, zpoles, zk) = _zpkdata(Hz, 1)

	slim = min(1e-3, max(1e-5, eps(1.0)^(1/(1+length(spoles)))))
	zlim = min(1e-3, max(1e-5, eps(1.0)^(1/(1+length(zzeros)))))

	H = zeros(size(f))
	w = 2*pi*f;	s = j*w;	z=exp.(s)
	for i in 1:length(f)
		wi = w[i];	si = s[i];	zi = z[i]
		cancel = false
		if !isempty(spoles)
			cancel = abs.(si .- spoles) .< slim
		end
		if !any(cancel)
			#wi is far from a pole, so just use the product Hs*Hz
			#H[i] = evalTF(Hs,si) * evalTF(Hz,zi);
			#use evalRPoly directly with z/p/k data, much faster. (David Alldred, Feb 4 2013)
			H[i] = sk * evalRPoly(szeros,si) / evalRPoly(spoles,si) *
				zk * evalRPoly(zzeros,zi) / evalRPoly(zpoles,zi)
		else
			#cancel pole(s) of Hs with corresponding zero(s) of Hz
			cancelz = abs.(zi .- zzeros) .< zlim
			sumcancel = sum(1*cancel); sumcancelz = sum(1*cancelz)
			if sumcancelz > sumcancel
				H[i] = 0.0
			elseif sumcancelz < sumcancel
				H[i] = Inf
			else
				H[i] = evalRPoly(szeros,si,sk) *
					zi^sumcancel * evalRPoly(zzeros[.!cancelz],zi,zk) /
					(evalRPoly(spoles[.!cancel],si,1)*evalRPoly(zpoles,zi,1))
			end
		end
	end
	return H
end
evalTFP(Hs,Hz,f) = evalTFP(Hs,Hz,[f])[1] #Single value


"""`evalMixedTF(tf,f, df=1e-5)`

Compute the mixed transfer function tf at a frequency f.
tf is a struct array with fields Hs and Hz wich represent 
a continuous-time/discrete-time tfs which must be multiplied together
and then added up.
"""
function evalMixedTF(tf::CTDTPrefilters, f; df=1e-5)
	H = zeros(size(f))
	for i = 1:length(tf)
		H = H + evalTFP(tf[i].Hs, tf[i].Hz, f)
	end

	err = findall(isnan.(H) .| isinf.(H))
	if !isempty(err)
		#Need to fill in the holes. !! Use a "nearby" frequency point
		for i in err
			H[i] = evalMixedTF(tf, f[i] .+ df, df=df*10)
		end
	end
	return H
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

"""`(pGain, Nimp) = powerGain(num,den,Nimp0=100)`

Calculate the power gain of a TF given in coefficient form.

 - `Nimp` is the recommended number of impulse response samples for use
   in future calls and `Nimp0` is the suggested number to use.
"""
function powerGain(num, den, Nimp0::Int=100)
	Nimp = Nimp0
	unstable = false

	sys = _tf(num,den,1)
	(imp,) = _impulse(sys, Nimp)
	if sum(abs.(imp[Nimp-10:Nimp])) < 1e-8 && Nimp > 50 #Use fewer samples.
		Nimp = round(Int, Nimp/1.3)
	else
		while sum(abs.(imp[Nimp-10:Nimp])) > 1e-6
			Nimp = Nimp*2
			(imp,) = _impulse(sys, Nimp)
			if sum(abs.(imp[Nimp-10:Nimp])) >= 50 || Nimp >= 1e4
				#H is close to being unstable
				unstable = true
				break
			end
		end
	end
	pGain = unstable ? Inf : sum(imp .^ 2)
	return (pGain, Nimp)
end

"""`zpk2 = cancelPZ(zpk1, tol=1e-6)`

Cancel zeros/poles in a zpk system
"""
function cancelPZ(zpk1, tol::Float64=1e-6)
	(z, p) = _zpkdata(zpk1, 1)
	zpk2 = deepcopy(zpk1) #Captures k, and whether tf is in s or z

	#Need to go in reverse order because z gets smaller with each cancellation
	for i = length(z):-1:1
		d = z[i] .- p
		cancel = findall(abs.(d) .< tol)
		if !isempty(cancel)
			deleteat!(p, cancel[1])
			deleteat!(z, i)
		end
	end
	zpk2.z = z
	zpk2.p = p
	return zpk2
end


#==infnorm
===============================================================================#
"""`nabsH(w, H)`

Computes the negative of the absolute value of H 
at the specified frequency on the unit circle. 

This function is used by infnorm().
"""
function nabsH(w, H)
	z = exp.(j*w)
	return -abs.(evalTF(H, z))
end


"""`(Hinf,fmax) = infnorm(H)`

Find the infinity norm of a z-domain transfer function.

Get a rough idea of the location of the maximum.
"""
function infnorm(H)
	N = 129
	w = range(0, stop=2π, length=N); dw = 2π/(N-1)
	Hval = abs.(evalTF(H, exp.(j*w)))
	wi = argmax(abs.(Hval))[1]
	Hinf = Hval[wi]

	local wmax
	if true #NEEDS_OPTIMIZATION
		opt = Optim.Options(x_tol=1e-8, f_tol=1e-6)
		f2min(w) = nabsH(w, H) #Function to minimize (optimum point)
		lower = w[wi]-dw; upper = w[wi]+dw
		wmax = optimize(f2min, lower, upper, GoldenSection()).minimum
	else
		msg = "Hinf: Warning. Optimization toolbox functions not found.\n" *
			" The result returned may not be very accurate."
		@warn(msg)
		wmax = w[wi]
	end

	Hinf = -nabsH(wmax, H)
	fmax = wmax/(2*pi)
	return (Hinf, wmax)
end


#==calculateTF
===============================================================================#
#Use a loophole to set complex poles and zeros in zpk objects
#MALaforge: I think this tries to implement minreal()
function setPolesAndZeros(z,p,k)
	tol = 3e-5 #tolerance for pole-zero cancellation
	z = deepcopy(z) #Don't change original
	Z = _zpk(ComplexF64[0.0],ComplexF64[],1.0,1.0)
	ztf = _zpk(ComplexF64[],ComplexF64[],k,1.0)
	for i in 1:length(p)
		match = abs.(p[i] .- z) .< tol
		if any(match)
			f = findall(match)
			deleteat!(z, f[1])
		else
			#ztf = ztf/(Z-p[i]) #Original code
			push!(ztf.p, p[i]) #VERIFYME: Think this is intent of original code
		end
	end
	for i in 1:length(z)
		#ztf = ztf*(Z-z[i]) #Original code
		push!(ztf.z, z[i]) #VERIFYME: Think this is intent of original code
	end
	return ztf
end

"""`(NTF,STF) = calculateTF(ABCD, k=1, wantSTF=true)`

Calculate the NTF and STF of a delta-sigma modulator whose loop filter
is described by the ABCD matrix, assuming a quantizer gain of k.
The NTF and STF are zpk objects.
"""
function calculateTF(ABCD, k::Vector=[1], wantSTF::Bool=true)
	nq = length(k)
	n = size(ABCD,1) - nq
	nu = size(ABCD,2) - size(ABCD,1)

	A,B,C,D = partitionABCD(ABCD, nu+nq)
	B1 = B[:,1:nu]
	B2 = B[:,nu+1:end]
	D1 = D[:,1:nu]
	if any(any( D[:,nu+1:end] .!= 0 ))
		throw("D2 must be zero")
	end
	K = diagm(k)
	#Find the noise transfer function by forming the closed-loop
	#system (sys_cl) in state-space form.
	Acl = A + B2*K*C
	Bcl = hcat(B1 + B2*K*D1, B2)
	Ccl = K*C
	Dcl = [K*D1 eye(nq)]
	#map(display, [Acl, Bcl, Ccl, Dcl])

	local NTF, STF
	tol = min(1e-3, max(1e-6, eps(1.0)^(1/(size(ABCD,1)))))
	if false # true || all(imag.(ABCD) .== 0) #real modulator
		#REQUESTHELP
		sys_cl = _ss(Acl,Bcl,Ccl,Dcl,1)
		sys_tf = ControlSystems.tf(sys_cl)
		mtfs = _minreal(sys_tf, tol)
		display(mtfs)
		@show (mtfs.nu, mtfs.ny), (nu, nq)
		#MALaforge: Can't make sense of the following; not sure what structure should be:
		STF = mtfs[:,1:nu]
		NTF = mtfs[:,nu+1:nu+nq]
	else #quadrature modulator and less advanced language features
		p = eig(Acl)
		NTFz = eig(A)
		NTF = setPolesAndZeros(NTFz, p, 1)
		if wantSTF
	        @warn("Sorry. calculateTF cannot compute the STF for a complex modulator with this version of the code.")
   	 end
	    STF=nothing
	end
	return (NTF, STF)
end


#Last line
