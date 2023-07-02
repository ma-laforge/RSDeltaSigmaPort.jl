#RSDeltaSigma: designLCBP() algorithm
#-------------------------------------------------------------------------------

#==Types
===============================================================================#
"""`LCBPParameters`

A struct containing the (`n`, `OSR`, `Hinf`, `f0`, `t` and `form`) arguments
for the LCBP filter, plus the following fields:
 - `t`: A two-element vector containing the start and end times of the
   feedback waveform. Default `= [0 1]`.
 - `l`: The inductances in each tank.
 - `c`: The capacitances in each tank.
 - `gu`: The transconductance from the u input to each of the tanks.
   The final term is the voltage gain from u to the comparator input.
   Default `= [1 0 .. ]`. (1 by n+1)
 - `gv`: The transconductance from the v output to each of the tanks. (1 by n+1)
 - `gw`: The inter-stage transconductances. Default `= [1 1 .. ]`. (1 by n-1)
 - `gx`: The gains from the top of each tank resistor to the comparator.
   Default `= [0 0 .. 1]`. (1 by n)
 - `rx`: The resistances inserted between the output of each interstage 
   transconductance and top of each tank. Default `= [0 ..]`. (1 by n)
   Note that rx[1] is not used.
 - `rl`: The series resistance of each inductor. Default `= [0 ..]`. (1 by n)
 - `gc`: The parallel conductances of each tank. Default `= [0 ..]`. (1 by n)
"""
mutable struct LCBPParameters
	n::Int
	OSR::Int
	Hinf::Float64
	f0::Float64
	t
	form::Symbol
	l::Array{Float64}; c::Array{Float64}
	gu::Array{Float64}; gv::Array{Float64}; gw::Array{Float64}; gx::Array{Float64}
	rx::Array{Float64}; rl::Array{Float64}; gc::Array{Float64}
end
LCBPParameters(n=3; OSR=64, Hinf=1.6, f0=0.25, t=[0 1], form::Symbol=:FB) =
	LCBPParameters(n, OSR, Hinf, f0, t, form, [], [], [], [], [], [], [], [], [])


#==
===============================================================================#
"""`(H,L0,ABCD,k)=LCparam2tf(param,k=1)`

Compute the transfer functions (and the Q. gain, k) of an LC modulator.

# Arguments
 - `param`: A LCBPParameters struct.
 - `k`: The value to assume for the quantizer gain.

# Note:
 - If `k` is omitted or if `k` is defaulted, `k` is taken to be 1.
 - If `k` is `0`, a value is computed by the formula
   `k = mean(abs(y))/mean(y.^2)`, where y is the quantizer input
   sequence when the modulator is fed with a small input at f0.

# Output
 - `H`: The closed-loop noise transfer function (in z).
 - `L0`: The open-loop tf (in s) from u to v.  G  may be evaluated
   using `evalContSTF(L0,H,f)`.
 - `k`: The value of the quantizer gain used to compute the above tfs.
 - `ABCD`: The state-space description of the discrete-time system

For the conversion from continuous to discrete time, corrections to the
standard MATLAB formulae need to be applied in order for the sampled-data
value to be taken at the END of the time interval.

This effectively makes t1>0 and requires corrections if t2>=1 and gv(n+1)~=0.
More complex formulae are needed for t2>=1, so for the moment this code
requires t2<1.
"""
function LCparam2tf(param::LCBPParameters, k::Float64=1.0)
	l = param.l; c = param.c
	f0 = mean(1 ./ sqrt.(l .* c)) ./ (2π *pi)
	n = length(l)

	gu = param.gu
		isempty(gu) && (gu = vcat(1, zeros(n))) #Defaults, if needed
	gu = padr(gu, n+1, 0) #Pad out with zeros, if needed.
	gv = padr(param.gv, n+1, 0) #Pad out with zeros, if needed.
	gw = param.gw
		isempty(gw) && (gw = ones(n-1)) #Defaults, if needed
	gx = param.gx
		isempty(gx) && (gx = vcat(zeros(n-1), 1)) #Defaults, if needed
	rx = param.rx
	if isempty(rx)
		rx = zeros(n)
	elseif 1==length(rx)
		throw("Not sure: Original code did not seem to make sense.")
		rx = ones(n)
	end
	t = param.t
		isempty(t) && (t = [0 1]) #Defaults, if needed
	rl = param.rl
	if isempty(rl)
		rl = zeros(n)
	elseif 1==length(rl)
		throw("Not sure: Original code did not seem to make sense.")
		rl = ones(n)
	end
	gc = param.gc
	if isempty(gc)
		gc = zeros(n)
	elseif 1==length(gc)
		throw("Not sure: Original code did not seem to make sense.")
		gc = ones(n)
	end

	#Form the state-space description of the modulator
	#The states are ordered v1, i1, v2, i2 ...
	n2 = 2*n
	local y
	ABc = zeros(n2,n2+2); CDc = zeros(1,n2+2)
	for i = 1:n
		i2 = 2*i-1
		#y represents the voltage at the top of the resistor stacked on the tank.
		if i == 1
			y = vcat(1, zeros(n2+1))
		else
			ABc[i2,:] = gw[i-1]*y
			y = rx[i]*gw[i-1]*y; y[i2] = 1
		end
		ABc[i2,n2+1:n2+2] = [gu[i] -gv[i]]
		ABc[i2:i2+1,i2:i2+1] = [-gc[i] -1; 1 -rl[i]]
		ABc[i2,:] = ABc[i2,:]/c[i]
		ABc[i2+1,:] = ABc[i2+1,:]/l[i]
		CDc = CDc + (gx[i]*y)'
	end
	CDc = CDc + [zeros(1,n2) gu[n+1] -gv[n+1]]

	#THIS IS THE NEW CODE !!
	(Ac,Bc,Cc,Dc) = partitionABCD([ABc;CDc], 2)
	sys_c = _ss( Ac, Bc[:,2], Cc, Dc[2] )
	(sys,) = mapCtoD(sys_c, t=t, f0=f0)
	#Augment the input matrix. 
	A=sys.A
	B=[padb(Bc[:,1],size(sys.B,1)) sys.B]
	C=sys.C
	D=[Dc[1] sys.D]

	#Compute L0; use the LC parameters to compute the poles and thereby
	#ensure exact agreement with the zeros in H.
	#!! THIS IS GOING TO NEED FIXING
	s_poles = zeros(2*n) .+ 0j #Make complex
	for i in 1:n
		_v = [l[i]*c[i], rl[i]*c[i] + gc[i]*l[i], 1+rl[i]*gc[i]]
		s_poles[2*i-1:2*i] = _roots(_v)
	end
	LF0 = _ss(Ac, Bc[:,1], Cc, Dc[1])
	L0 = _zpk( LF0 )
	L0.p = s_poles

	#Compute H. Use the poles in L0 to compute the zeros in H.
	ABCD =[A B; C D]
	if k==0 #Estimate k by simulatiing the modulator with a small input
		w0 = mean(1/sqrt.(l .* c))
		H = calculateTF(ABCD)
		H.z = exp.(s_poles)
		stf0 = abs(evalTFP(L0, H, w0/(2π)))
		u = 0.1/stf0*sin.(w0*(0:10000))
		y = simulateDSM(u, [A B; C D]).y
		k = mean(abs.(y))/mean(y .^ 2)
	end
	(H, STF) = calculateTF(ABCD, [k])
	H.z = exp.(s_poles)

	#Correct L0k to include the quantizer gain
	L0.k = L0.k*k
	return (H,L0,ABCD,k)
end

#Dead code???
function XLCoptparam2tf(x, param)
	n = param.n
	#Uncomment the following line to take into account the quantizer gain 
	#before returning results and doing the plots
	#(H, L0, ABCD, k) = LCparam2tf(param, 0)
	#Uncomment the following lines to yield parameters corresponding to k=1
	#param.gu = k*param.gu; param.gv = k*param.gv
	#ABCD[2*n+1, :] = ABCD[2*n+1, :] / k
	#ABCD[:, 2*n+[1 2]] = ABCD[:, 2*n+[1 2]] * k

	#Now fix up the inputs for 0dB gain at f0.
	gain_f0 = abs(evalTFP(L0, H, f0))
	param.gu = param.gu/gain_f0; L0.k = L0.k/gain_f0
	ABCD[:,2*n+1] = ABCD[:,2*n+1]/gain_f0
	if (:FF==form) || (:FFR==form)
		param.gu[n+1] = 1 #For the flattest STF
	end
	if dbg
		@warn("LCplotTF: not yet implemented")
		#LCplotTF(H,L0,param)
	end
	return (H, L0, ABCD, param)
end

"""`LCoptparam2tf(...)`

Convert optimization parameters to transfer functions.
"""
function LCoptparam2tf(x,param)
	n = param.n
	param.gw = ones(n-1)

	if :FB == param.form
		param.gu = vcat(1, zeros(n))
		param.gv = vcat(x)
		param.gx = vcat(zeros(n-1), 1)
		param.rx = zeros(n)
	elseif :FF == param.form
		param.gu = vcat(1, zeros(n-1), 1)
		param.gv = vcat(1, zeros(n))
		param.gx = vcat(x)
		param.rx = zeros(n)
	elseif :FFR == param.form
		param.gu = vcat(1, zeros(n-1), 1)
		param.gv = vcat(x[1], zeros(n))
		param.gx = vcat(zeros(n-1), 1)
		param.rx = vcat(0, x[2:n])
	elseif :GEN == param.form
		param.gv = x[1:n]
		param.rx = vcat(0, x[n+1:2*n-1])
	else
		form = param.form
		throw("form=:$form is not supported.")
	end

	(H, L0, ABCD) = LCparam2tf(param)
	return (H,L0,ABCD,param)
end


"""`f = LCObj1a(x, param, max_radius, dbg)`

The objective function for the initial optimization process used to put
the roots of the denominator of the LCBP NTF inside the unit circle.
"""
function LCObj1a(x, param, max_radius, dbg)
	return f=1 #No actual objective
end

"""`(C, Ceq) = LCObj1b(x, param, max_radius, dbg)`

The constraints function for the initial optimization process used to put
the roots of the denominator of the LCBP NTF inside the unit circle.
"""
function LCObj1b(x, param, max_radius, dbg)
	(H,) = LCoptparam2tf(x,param)
	objective = 1 #No actual objective

	rmax = maximum(abs.(H.p[:]))
	C = rmax - max_radius
	Ceq = []
	if dbg
		@show(x)
		@show(rmax)
	end
	return (C, Ceq)
end

#The objective function
function LCObj2a(x, param, rmax, dbg)
	if dbg
		@show(x)
	end
	(H, L0, ABCD) = LCoptparam2tf(x, param)
	#The objective is the in-band noise, in dB
	f1 = param.f0 - 0.25/param.OSR
	f2 = param.f0 + 0.25/param.OSR
	objective = dbv(rmsGain(H,f1,f2)/sqrt(3*param.OSR))
	println(@sprintf("N0^2 = %.1fdB", objective))
	return objective
end

"""`(C, Ceq) = LCObj2b(x,param,rmax,dbg)`
% The constraints: r-rmax, Hinf-Hinfdesired, quantizer input,
% and for the FB form, |stf|-1.1.
"""
function LCObj2b(x, param, rmax, dbg)
	(H, L0, ABCD) = LCoptparam2tf(x, param)
	Ceq = []
	local C

	#The constraints are on the maximum radius of the poles of the NTF
	#the infinity-norm of the NTF, the peaking of the STF (for the 'FB' form),
	#and the maximum input to the quantizer.
	max_radius = maximum(abs(H.p[:]))
	H_inf = infnorm(H)
	stf0 = abs( evalTFP(L0, H, param.f0) )
	y = simulateDSM(0.5/stf0*sin(2π*param.f0*(0:1000)), ABCD).y
	ymax = maximum(abs.(y))
	if :FB == param.form
		f1 = param.f0 - 0.25/param.OSR
		f2 = param.f0 + 0.25/param.OSR
		stf1 = abs( evalTFP(L0, H, f1) )/stf0
		stf2 = abs( evalTFP(L0, H, f2) )/stf0
		C = vcat(100*(max_radius .- rmax), H_inf-param.Hinf, stf1-1.01, stf2-1.01, (ymax-5)/10)
		#C = vcat(100*(max_radius .- rmax), H_inf-param.Hinf, stf1-1.01, stf2-1.01, log(ymax/5))
	else
		C = vcat(100*(max_radius .- rmax), H_inf-param.Hinf, (ymax-5)/10)
	end

	if dbg
		println(@sprintf("constraints = %8.5f", C))
				println(@sprintf("rmax = %5.3f, Hinf = %4.2f, ymax = %5.3f\n",
					max_radius, H_inf, ymax)
				)
	end
	return (C, Ceq)
end


#==Main algorithm
===============================================================================#
"""`(param,H,L0,ABCD,x)=designLCBP(n=3,OSR=64,opt=2,Hinf=1.6,f0=1/4,t=[0 1],form='FB',x0,dbg=0)`

Design an LC modulator of order `2*n` with center frequency `f0`.
The feedback waveform is a rectangular pulse from `t[1]` to `t[2]`.

# Arguments
 - `n`: The number of LC tanks.
 - `OSR`: The oversampling ratio.
 - `opt`: A flag for NTF zero optimization. See synthesizeNTF.
 - `Hinf`: The target out-of-band gain of the NTF.
 - `f0`: The center frequency of the modulator, normalized to fs.
 - `t`: A 2-element vector containing the (normalized) feedback pulse edge times.
 - `form`: The modulator topology. One of `:FB`, `:FF`, `:FFR` or `:GEN`.
 - `x0`: An initial guess at the parameter vector for the optimizer.
 - `dbg`: A flag for enabling debugging/progress report output.

# Output
 - `param`: A LCBPParameters struct.
 - `H`: The noise transfer function of the modulator, taking into account
   the post-design quantizer gain.
 - `L0`: A description of the open-loop TF from u to v, for use with evalTFP().
 - `ABCD`: The state-space representation of the discrete-time equivalent.
"""
function designLCBP(n=3; OSR=64, opt=2, Hinf=1.6, f0=0.25, t=[0 1],
		form::Symbol=:FB, x0=NaN, dbg::Bool=false)
	param = LCBPParameters(n, OSR=OSR, Hinf=Hinf, f0=f0, t=t, form=form)
	#Compute the resonant frequencies for the LC tanks 
	#and the L and C values.
	local w
	if opt>0
		w = 2π*(f0 .+ 0.25/OSR*ds_optzeros(n,opt))
	else
		w = 2π*f0*ones(n)
	end
	param.l = (1 ./ w); param.c = (1 ./ w) #Impedance scaling is possible.
	#Find values for the parameters that yield a stable system.
	#?? Could some control theory help here ??
	x = x0
	if isnan(x0)
		x = collect(2:n+1)*1.0 #Needs to be Float64
	elseif length(x0) != n
		throw("Initial guess (x0) has the wrong number of parameters.")
	end
	if dbg
		println("Iterations for initial pole placement:")
	end
	max_radius = 0.97
	opt = Optim.Options(x_tol=0.001, f_tol=1, iterations=1000)
	lx = Float64[0,0,0]; ux = Float64[10,10,10]
	lc = [-Inf]; uc = [max_radius^2]
	fun1!(x) = LCObj1a(x, param, max_radius, dbg)
	con_c1!(c, x) = ((C, Ceq)=LCObj1b(x, param, max_radius, dbg); C)
#	dfc = TwiceDifferentiableConstraints(con_c!, con_jacobian!, con_h!,lx,ux,lc,uc)
	dfc = TwiceDifferentiableConstraints(con_c1!,lx,ux,lc,uc)#,:forward)
#	x = optimize(fun1!, dfc, x, IPNewton()).minimum
	x = optimize(fun1!, con_c1!, x, LBFGS()).minimum
	#x = fmincon(@LCObj1a,x,[],[],[],[],[],[],@LCObj1b,options,param,0.97,dbg);
	(H,) = LCoptparam2tf(x,param)
	rmax = maximum(abs.(H.p[:]))
	if rmax>max_radius #Failure to converge!
		throw("Unable to find parameter values which stabilize the system.")
	end
	#Do the main optimization for the feedback parameters.
	if dbg
		println("\nParameter Optimization:")
	end
	max_radius = 0.97
	opt = Optim.Options(x_tol=0.001, f_tol=1, iterations=1000,
		show_trace=dbg, #Extra debugging info
		g_tol=0.01 #Solver-specific
	)
	lx = Float64[]; ux = Float64[]
	lc = [-Inf]; uc = [max_radius^2]
	fun2!(x) = LCObj2a(x, param, max_radius, dbg)
	con_c2!(c, x) = ((C, Ceq)=LCObj2b(x, param, max_radius, dbg); C)
	dfc = TwiceDifferentiableConstraints(con_c2!,lx,ux,lc,uc,:forward)
	x = optimize(fun2!, dfc, x, IPNewton()).minimum
	#x = fmincon(@LCObj2a,x,[],[],[],[],[],[],@LCObj2b,options,param,0.97,dbg);
	return (param,H,L0,ABCD,x)
end

#Last line
