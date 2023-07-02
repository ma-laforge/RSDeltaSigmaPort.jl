#RSDeltaSigmaPort: Synthesize NTFs
#-------------------------------------------------------------------------------


#==Warnings
===============================================================================#
function _synthesizeNTF_warn_Hinf()
	msg = string(@__FILE__, ": Unable to achieve specified Hinf.\n",
		"Setting all NTF poles to zero."
	)
	@warn(msg)
end
function _synthesizeNTF_warn_iter()
	msg = string(@__FILE__, ": Danger! Iteration limit exceeded.")
	@error(msg)
end
function _synthesizeNTF_warn_Hinf_iter()
	msg = string(@__FILE__, ": Danger! Hinf iteration limit exceeded.")
	@error(msg)
end


#==synthesizeNTF_setup()
===============================================================================#
function synthesizeNTF_setup(order::Int, OSR::Int, opt, f0::Float64, isNTF1::Bool)
	local z
	order = Float64(order) #Type stability
	dw = π/OSR

	if f0!=0		#Bandpass design-- halve the order temporarily.
		order = order/2
		dw = π/(2*OSR)
	end
	order_i = array_round(order)

	#Determine the zeros.
	if length(opt)==1
		if opt==0
			z = zeros(order_i)
		else
			opt_for_zeros = isNTF1 ? (1+mod(opt-1,2)) : opt
			z = dw*ds_optzeros(order, opt_for_zeros)
			if isempty(z)
				throw_unreachable()
				return;
			end
		end
		if f0!=0		# Bandpass design-- shift and replicate the zeros.
			order = order*2
			z = z .+ 2π*f0
			z = conj.(z)
			z = vcat(z, -z)
		end
		z = exp.(j*z)
	else
		z = opt[:] * complex(1.0) #Ensure is Vector{Complex{Float64}}
	end

	#Order might possibly be a N.5 value (float) here. Not sure if this is actually the intent.
	return (order, z) #Order potentially gets modified. Not sure this is actually the intent.
end


#==synthesizeNTF0()
===============================================================================#
function _synthesizeNTF0_z(order::Float64, OSR::Int, z::Vector{Complex{Float64}}, H_inf::Float64, f0::Float64)
	order_i = array_round(order)
	Ts = 1.0
	ntf = _zpk(z,zeros(order_i),1,Ts)
	itn_limit = 100
	#ntf_z, ntf_p, ntf_k, = zpkdata(ntf)

	# Iteratively determine the poles by finding the value of the x-parameter
	# which results in the desired H_inf.
	if f0 == 0			# Lowpass design
		HinfLimit = 2^order  # !!! The limit is actually lower for opt=1 and low OSR
		if H_inf >= HinfLimit
			_synthesizeNTF_warn_Hinf()
			ntf.p = zeros(order_i)
		else
			x=0.3^(order-1) # starting guess
			converged = 0
			local fprev = 0 #Initialize
			local delta_x = 0 #Initialize
			for itn=1:itn_limit
				me2 = -0.5*(x^(2/order))
				w = (2*collect(1:order_i) .- 1)*π/order #VERIFYME
				mb2 = 1 .+ me2*exp.(j*w)
				p = mb2 - sqrt.(mb2 .^ 2 .- 1)
				out = findall(abs.(p) .> 1)
				p[out] = 1 ./ p[out] # reflect poles to be inside the unit circle.
				ntf.z = z; ntf.p = cplxpair(p)
				f = real(evalTF(ntf,-1))-H_inf
				# [ x f ]
				if itn==1 
					delta_x = -f/100
				else
					delta_x = -f*delta_x/(f-fprev)
				end

				xplus = x+delta_x
				if xplus>0 
					x = xplus
				else
					x = x*0.1
				end
				fprev = f;

				if abs(f)<1e-10 || abs(delta_x)<1e-10 
					converged = 1
					break
				end
				if x>1e6
					_synthesizeNTF_warn_Hinf()
					ntf.z = z; ntf.p = zeros(order_i)
					break;
				end
				if itn == itn_limit
					_synthesizeNTF_warn_iter()
				end
			end
		end
	else # Bandpass design.
		x = 0.3^(order/2-1)	# starting guess (not very good for f0~0)
		if f0>0.25
			z_inf=1
		else
			z_inf=-1
		end
		c2pif0 = cos(2π*f0)
		local fprev = 0 #Initialize
		local delta_x = 0 #Initialize
		for itn=1:itn_limit
			e2 = 0.5*x^(2/order)
			w = (2*collect(1:order_i) .- 1)*π/order #VERIFYME

			mb2 = c2pif0 .+ e2*exp.(j*w)
			p = mb2 .- sqrt.(mb2.^2 .- 1)
			# reflect poles to be inside the unit circle.
			out = findall(abs.(p) .> 1)
			p[out] = 1 ./ p[out]
			ntf.z = z; ntf.p = cplxpair(p)
			f = real(evalTF(ntf,z_inf))-H_inf
			# [x f]
			if itn==1 
				delta_x = -f/100
			else
				delta_x = -f*delta_x/(f-fprev)
			end
	
			xplus = x+delta_x
			if xplus > 0
				x = xplus
			else
				x = x*0.1
			end

			fprev = f
			if abs(f)<1e-10 || abs(delta_x)<1e-10
				break
			end
			if x>1e6
				_synthesizeNTF_warn_Hinf()
				ntf.p = zeros(order_i)
				break
			end
			if itn == itn_limit
				_synthesizeNTF_warn_iter()
			end
		end
	end

	ntf.z = cleancplxpair(ntf.z)
	ntf.p = cleancplxpair(ntf.p)
	return ntf
end

function synthesizeNTF0(order::Int, OSR::Int, opt, H_inf::Float64, f0::Float64)
	(order, z) = synthesizeNTF_setup(order, OSR, opt, f0, false)
	_synthesizeNTF0_z(order, OSR, z, H_inf, f0)
end


#==synthesizeNTF1()
===============================================================================#
function _synthesizeNTF1_z(order::Float64, osr::Int, opt, z::Vector{Complex{Float64}}, H_inf::Float64, f0::Float64)
	order_i = array_round(order)
	zp = z[findall(angle.(z) .> 0)]
	x0 = (angle.(zp) .- 2π*f0) * osr / π
	if opt==4 && f0 != 0
		#Do not optimize the zeros at f0
		deleteat!(x0, abs.(x0) .< 1e-10)
	end

	Ts = 1.0
	ntf = _zpk(z,zeros(order_i),1,Ts)
	Hinf_itn_limit = 100

	opt_iteration = 5 #Max number of zero-optimizing/Hinf iterations
	while opt_iteration > 0
		#Iteratively determine the poles by finding the value of the x-parameter
		#which results in the desired H_inf.
		ftol = 1e-10
		z_inf = f0>0.25 ? 1 : -1
		if f0 == 0 #Lowpass design
			HinfLimit = 2^order #!!! The limit is actually lower for opt=1 and low osr
			if H_inf >= HinfLimit
				_synthesizeNTF_warn_Hinf()
				ntf.p = zeros(order_i,1)
			else
				x=0.3^(order-1) #starting guess
				converged = false
				local delta_x
				local fprev
				for itn in 1:Hinf_itn_limit
					me2 = -0.5*(x^(2/order))
					w = (2*collect(1:order_i) .- 1)*π/order #VERIFYME
					mb2 = 1 .+ me2*exp.(j*w)
					p = mb2 .- sqrt.(mb2 .^ 2 .- 1)
					out = findall(abs.(p) .> 1)
					p[out] = 1 ./ p[out] #reflect poles to be inside the unit circle.
					p = cplxpair(p)
					ntf.z = z;	ntf.p = p
					f = real(evalTF(ntf, z_inf))-H_inf
					#[ x f ]
					if itn==1
						delta_x = -f/100
					else
						delta_x = -f*delta_x/(f-fprev)
					end
					xplus = x+delta_x
					x = (xplus>0) ? xplus : x*0.1
					fprev = f
					if abs(f)<ftol || abs(delta_x)<1e-10
						converged = true
						break
					end
					if x>1e6
						_synthesizeNTF_warn_Hinf()
						ntf.z = z;	ntf.p = zeros(order_i,1)
						break
					end
					if itn == Hinf_itn_limit
						_synthesizeNTF_warn_iter()
					end
				end
			end
		else #Bandpass design.
			x = 0.3^(order/2-1) #starting guess (not very good for f0~0)
			c2pif0 = cos(2π*f0)
			local delta_x
			local fprev
			for itn=1:Hinf_itn_limit
				e2 = 0.5*x^(2/order)
				w = (2*collect(1:order_i) .- 1)*π/order
				mb2 = c2pif0 .+ e2*exp.(j*w)
				p = mb2 .- sqrt.(mb2 .^ 2 .- 1)
				#reflect poles to be inside the unit circle.
				out = findall(abs.(p) .> 1)
				p[out] = 1 ./ p[out]
				p = cplxpair(p)
				ntf.z = z; ntf.p = p
				f = real(evalTF(ntf, z_inf)) - H_inf
				#[x f]
				if itn==1
					delta_x = -f/100
				else
					delta_x = -f*delta_x/(f-fprev)
				end
				xplus = x+delta_x
				x = (xplus>0) ? xplus : x*0.1
				fprev = f
				if abs(f)<ftol || abs(delta_x)<1e-10
					break
				end
				if x>1e6
					_synthesizeNTF_warn_Hinf()
					p = zeros(order_i,1)
					ntf.p = p
					break
				end
				if itn == Hinf_itn_limit
					_synthesizeNTF_warn_Hinf_iter()
				end
			end
		end
		########################################################
		if opt < 3 #Do not optimize the zeros
			opt_iteration = 0
		else
			if f0 == 0
				ub = ones(size(x0))
				lb = zeros(size(x0))
			else
				ub = 0.5*ones(size(x0))
				lb = -ub
			end
			opt = Optim.Options(x_tol=0.001, f_tol=0.01, iterations=100)
			f2min(x) = ds_synNTFobj1(x,p,osr,f0) #Function to minimize (optimum point)
			x = optimize(f2min, x0, lb, ub, Fminbox(GradientDescent())).minimum
			x0 = x
			z = exp.(2π*j*(f0 .+ 0.5/osr*x))
			if f0>0
				z = padt(z, array_round(length(p)/2), val=exp(2π*j*f0))
			end
			z = vcat(z, conj.(z))
			if f0==0
				z = padt(z, length(p), val=1.0)
			end
			ntf.z = z; ntf.p = p
			if abs( real(evalTF(ntf, z_inf)) - H_inf ) < ftol
				opt_iteration = 0
			else
				opt_iteration = opt_iteration - 1
			end
		end
	end

	ntf.z = cleancplxpair(ntf.z)
	ntf.p = cleancplxpair(ntf.p)
	return ntf
end

function synthesizeNTF1(order::Int, OSR::Int, opt, H_inf::Float64, f0::Float64)
	(order, z) = synthesizeNTF_setup(order, OSR, opt, f0, true)
	_synthesizeNTF1_z(order, OSR, opt, z, H_inf, f0)
end


#==synthesizeNTF()
===============================================================================#
"""`NTF = synthesizeNTF(order=3, OSR=64, opt=0, H_inf=1.5, f0=0)`

Synthesize a noise transfer function for a delta-sigma modulator.

# Inputs
 - `order`: order of the modulator
 - `OSR`: oversampling ratio
 - `opt`: flag for optimized zeros\r
 - `--> opt=0` -> not optimized,
 - `--> opt=1` -> optimized,
 - `--> opt=2` -> optimized with at least one zero at band-center,
 - `--> opt=3` -> optimized zeros (Requires MATLAB6 and Optimization Toolbox),
 - `--> opt=[z]` -> zero locations in complex form (TODO?)
 - `H_inf`: maximum NTF gain
 - `f0`: center frequency (1->fs)

# Output
`NTF` is a `_zpk` object containing the zeros and poles of the NTF
(See [`_zpk`](@ref)).

# Note

This is actually a wrapper function which calls either the 
appropriate version of `synthesizeNTF` based on the availability
of the `fmincon` function from the Optimization Toolbox
"""
function synthesizeNTF(order::Int=3, osr::Int=64; opt=0, H_inf::Float64=1.5, f0::Float64=0.0)
	if f0 > 0.5
		throw(ErrorException("f0 must be less than 0.5."))
	end
	if f0 != 0 && f0 < 0.25/osr
		@warn(string("(", @__FILE__, ") Creating a lowpass ntf."))
		f0 = 0
	end
	if f0 != 0 && mod(order,2) != 0
		throw(ErrorException("order must be even for a bandpass modulator."))
	end

	if length(opt)>1 && length(opt)!=order
		throw(ErrorException("The opt vector must be of length $order(=order)."))
	end

	if true #Use version that applies optimization
		#throw_notimplemented()
		return synthesizeNTF1(order,osr,opt,H_inf,f0)
	else
		return synthesizeNTF0(order,osr,opt,H_inf,f0)
	end
end

"""`NTF = synthesizeNTF(dsm)`

Synthesize a noise transfer function for a ΔΣ modulator.

# Output
`NTF` is a `_zpk` object containing the zeros and poles of the NTF
(See [`_zpk`](@ref)).

# See also:
 - `dsm` types: [`RealDSM`](@ref), [`QuadratureDSM`](@ref)
 - [`clans()`](@ref): "Closed-loop analysis of noise-shaper." An alternative
   method for selecting NTFs based on the 1-norm of the impulse
   response of the NTF.
 - [`synthesizeChebyshevNTF`](@ref): Select a type-2 highpass Chebyshev NTF.
   This function does a better job than `synthesizeNTF` if `OSR` or `H_inf`
   is low.
 - [`_zpk`](@ref)
""" synthesizeNTF

synthesizeNTF(dsm::RealDSM) =
	synthesizeNTF(dsm.order, dsm.OSR, opt=dsm.opt, H_inf=dsm.Hinf, f0=dsm.f0)

#synthesizeNTF(dsm::QuadratureDSM) =
#	synthesizeQNTF(dsm.order, dsm.OSR, dsm.f0, dsm.NG, dsm.ING)

#Last line
