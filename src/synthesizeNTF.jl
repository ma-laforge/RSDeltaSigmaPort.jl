#RSDeltaSigmaPort: Synthesize NTFs
#-------------------------------------------------------------------------------

function synthesizeNTF0_z(order::Float64, OSR::Int, z::Vector{Complex{Float64}}, H_inf::Float64, f0::Float64)
	order_i = array_round(order)
	Ts = 1.0
	ntf = _zpk(z,zeros(order_i),1,Ts)
	itn_limit = 100
	#ntf_z, ntf_p, ntf_k, = zpkdata(ntf)

	function warn_Hinf()
		msg = string(@__FILE__, ": Unable to achieve specified Hinf.\n",
			"Setting all NTF poles to zero."
		)
		@warn(msg)
	end
	function warn_iter()
		msg = string(@__FILE__, ": Danger! Iteration limit exceeded.")
		@error(msg)
	end

	# Iteratively determine the poles by finding the value of the x-parameter
	# which results in the desired H_inf.
	if f0 == 0			# Lowpass design
		HinfLimit = 2^order  # !!! The limit is actually lower for opt=1 and low OSR
		if H_inf >= HinfLimit
			warn_Hinf()
			ntf.p = zeros(order_i)
		else
			x=0.3^(order-1) # starting guess
			converged = 0
			local fprev = 0 #Initialize
			local delta_x = 0 #Initialize
			for itn=1:itn_limit
				me2 = -0.5*(x^(2/order))
				w = (2*collect(1:order_i) .- 1)*pi/order #VERIFYME
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
					warn_Hinf()
					ntf.z = z; ntf.p = zeros(order_i)
					break;
				end
				if itn == itn_limit
					warn_iter()
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
			w = (2*collect(1:order_i) .- 1)*pi/order #VERIFYME

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
				warn_Hinf()
				ntf.p = zeros(order_i)
				break
			end
			if itn == itn_limit
				warn_iter()
			end
		end
	end

	ntf.z = cleancplxpair(ntf.z)
	ntf.p = cleancplxpair(ntf.p)
	return ntf
end

function synthesizeNTF0(order::Int, OSR::Int, opt, H_inf::Float64, f0::Float64)
	local z
	order = Float64(order) #Type stability
	dw = pi/OSR

	if f0!=0		#Bandpass design-- halve the order temporarily.
		order = order/2
		dw = pi/(2*OSR)
	end
	order_i = array_round(order)

	#Determine the zeros.
	if length(opt)==1
		if opt==0
			z = zeros(order_i)
		else
			z = dw*ds_optzeros(order,opt)
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

	return synthesizeNTF0_z(order, OSR, z, H_inf, f0)
end


"""`ntf = synthesizeNTF(order=3,osr=64,opt=0,H_inf=1.5,f0=0)`
Synthesize a noise transfer function for a delta-sigma modulator.
	order =	order of the modulator
	osr =	oversampling ratio
	opt =	flag for optimized zeros
		0 -> not optimized,
		1 -> optimized, 
		2 -> optimized with at least one zero at band-center
		3 -> optimized zeros (Requires MATLAB6 and Optimization Toolbox)
       [z] -> zero locations in complex form
	H_inf =	maximum NTF gain
	f0 =	center frequency (1->fs)

ntf is a zpk object containing the zeros and poles of the NTF. See zpk.m

 See also 
  clans()   "Closed-loop analysis of noise-shaper." An alternative
            method for selecting NTFs based on the 1-norm of the 
            impulse response of the NTF

  synthesizeChebyshevNTF()    Select a type-2 highpass Chebyshev NTF.
            This function does a better job than synthesizeNTF if osr
            or H_inf is low.

 This is actually a wrapper function which calls either the 
 appropriate version of synthesizeNTF based on the availability
 of the 'fmincon' function from the Optimization Toolbox
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

	if FMINCON_AVAIL
		throw_notimplemented()
		#return synthesizeNTF1(order,osr,opt,H_inf,f0)
	else
		return synthesizeNTF0(order,osr,opt,H_inf,f0)
	end
end

synthesizeNTF(dsm::RealDSM) =
	synthesizeNTF(dsm.order, dsm.OSR, opt=dsm.opt, H_inf=dsm.Hinf, f0=dsm.f0)

#synthesizeNTF(dsm::QuadratureDSM) =
#	synthesizeQNTF(dsm.order, dsm.OSR, dsm.f0, dsm.NG, dsm.ING)

#Last line
