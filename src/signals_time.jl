#RSDeltaSigma: time-domain signal generators
#-------------------------------------------------------------------------------


#==Input signal generators
===============================================================================#
"""`genTestTone_sin()`

Limiation:
If N is small relative to ftest (you don't generate a whole period), iftest gets rounded to nothing.
Typical ampdB: -3dB
"""
function genTestTone_sin(ampdB::Float64, ftest::Float64; M::Int=1, phi0::Float64=0.0, N::Int=2^12)
	amp = undbv(ampdB) #Test tone amplitude, relative to full-scale.
	iftest = round(Int, ftest*N)
	u = amp*M*sin.((2π/N)*iftest * (0:N-1) .+ phi0)
	return (u, iftest)
end

function genTestTone_quad(ampdB::Float64, ftest::Float64; M::Int=1, phi0::Float64=0.0, N::Int=2^12)
	amp = undbv(ampdB) #Test tone amplitude, relative to full-scale.
	iftest = round(Int, ftest*N)
	u = amp*M*exp.((2π*j/N)*iftest * (0:N-1) .+ phi0)
	return (u, iftest)
end

genTestTone(dsm::RealDSM, ampdB::Float64, ftest::Float64; phi0::Float64=0.0, N::Int=2^12) =
	genTestTone_sin(ampdB, ftest, M=dsm.M, phi0=phi0, N=N)
genTestTone(dsm::QuadratureDSM, ampdB::Float64, ftest::Float64; phi0::Float64=0.0, N::Int=2^12) =
	genTestTone_quad(ampdB, ftest, M=dsm.M, phi0=phi0, N=N)


#==impL1
===============================================================================#
"""`y=impL1(ntf, n=10)`

Compute the impulse response from the comparator
output to the comparator input for the given NTF.
n is the (optional) number of points (10).

This function is useful when verifying the realization
of a NTF with a specified topology.
"""
function impL1(ntf::ZPKData, n::Int=10)
	(z, p) = _zpkdata(ntf)

	lf_den = padr(poly(z), length(p)+1)
	lf_num = lf_den .- poly(p)
	if any(imag.(vcat(lf_num, lf_den)) .!= 0)
		#Complex loop filter
		lfr_den = real.( conv(lf_den, conj.(lf_den)) )
		lfr_num = conv(lf_num, conj.(lf_den))
		lf_i = _tf(real.(lfr_num), lfr_den, 1)
		lf_q = _tf(imag.(lfr_num), lfr_den, 1)
		(yi,) = _impulse(lf_i, n)
		(yq,) = _impulse(lf_q, n)
		y = yi .+ j*yq
	else
		lf_num = real.(lf_num); lf_den = real(lf_den)
		(y,) = _impulse(_tf(lf_num, lf_den, 1), n)
	end
	return y
end


#==pulse
===============================================================================#
"""`y = pulse(S, tp=[0 1], dt=1, tfinal=10, nosum=true)`

Calculate the sampled pulse response of a ct system. `tp` may be an array of
pulse timings, one for each input.

# Outputs
 - `y` The pulse response.

# Inputs
 - `S`: An LTI object specifying the system.
 - `tp`: An nx2 array of pulse timings.
 - `dt`: The time increment.
 - `tfinal`: The time of the last desired sample.
 - `nosum`: A flag indicating that the responses are not to be summed.
"""
function pulse(S, tp=[0 1], dt::Float64=1.0, tfinal::Float64=10.0, nosum::Bool=false)
	#Check the arguments
	if isdiscrete(S)
		error("S must be a cts-time system.")
	end

	#Compute the time increment
	dd = 1
	for i in 1:prod(size(tp))
		(x, di) = rat(tp[i],1e-3)
		dd = lcm(di,dd);
	end
	(x, ddt) = rat(dt,1e-3)
	(x, df) = rat(tfinal,1e-3)
	delta_t = 1 / lcm( dd, lcm(ddt,df) )
	delta_t = max(1e-3, delta_t) #Put a lower limit on delta_t
	(y1, ty1, xy1) = step(S, 0:delta_t:tfinal)

#= #DEBUG
	p=cons(:plot, title="")
	for i in 1:size(y1, 3)
	push!(p, cons(:wfrm, waveform(collect(ty1[:]), y1[:,1,i])))
	end
	displaygui(p)
=#

	nd = round(Int, dt/delta_t)
	nf = round(Int, tfinal/delta_t)
	ndac = size(tp,1)
	ni = size(S.B,2)
	if mod(ni,ndac) != 0 
		error("The number of inputs must be divisible by the number of dac timings.")
		#This requirement comes from the complex case, where the number of inputs
		#is 2 times the number of dac timings. I think this could be tidied up.
	end
	nis = div(ni,ndac) #Number of inputs grouped together with a common DAC timing
	# (2 for the complex case)
	if !nosum #Sum the responses due to each input set
		y = zeros(array_round(tfinal/dt+1),size(S.C,1),nis)
	else
		y = zeros(array_round(tfinal/dt+1),size(S.C,1),ni)
	end

	for i = 1:ndac
		n1 = round(Int, tp[i,1]/delta_t)
		n2 = round(Int, tp[i,2]/delta_t)
		z1 = (n1, size(y1,2), nis)
		z2 = (n2, size(y1,2), nis)
		yy = [zeros(z1); y1[1:nf-n1+1,:,(i-1)*nis+1:i*nis]] .-
		     [zeros(z2); y1[1:nf-n2+1,:,(i-1)*nis+1:i*nis]]
		yy = yy[1:nd:end,:,:]
		if !nosum #Sum the responses due to each input set
			y = y .+ yy
		else
			y[:,:,i] = yy
		end
	end

	return y
end


#==pulse
===============================================================================#
"""`sv = ds_therm(v,M)`

`sv` is an `M` by `length(v)` matrix of `+/-1`s
where `sum(sv)=v` and the `-1`s are located in the bottom rows
"""
function ds_therm(v, M::Int)
	sv = -ones(M, length(v))
	for i in 1:length(v)
		iend = array_round((M+v[i])/2)
		sv[1:iend, i] .= 1
	end
	return sv
end

#Last line
