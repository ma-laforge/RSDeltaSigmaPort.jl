#RSDeltaSigma: realize continuous-time NTF
#-------------------------------------------------------------------------------

#==Helper functions
===============================================================================#


#==realizeNTF_ct
===============================================================================#
"""`(ABCDc, tdac2) = realizeNTF_ct(NTF, form='FB', tdac, ordering=[1:n], bp=zeros(...), ABCDc)`

Realize a NTF with a continuous-time loop filter.

# Output
 - `ABCDc`: A state-space description of the CT loop filter
 - `tdac2`: A matrix with the DAC timings, including ones that were
   automatically added.

# Input Arguments
 - `NTF`: A noise transfer function in pole-zero form.
 - `form`: `= {:FB, :FF}` Specifies the topology of the loop filter.
 - `-->` For the `:FB` structure, the elements of `Bc` are calculated
   so that the sampled pulse response matches the `L1` impulse respnse.
 - `-->` For the `:FF` structure, `Cc` is calculated.
 - `tdac`: The timing for the feedback DAC(s).
 - `-->` If `tdac[1]>=1`, direct feedback terms are added to the quantizer.
 - `-->` Multiple timings (1 or more per integrator) for the `:FB`
   topology can be specified by making `tdac` a cell array, e.g.
   `tdac = { [1,2]; [1 2]; {[0.5 1],[1 1.5]}; []}`.
 - `-->` In above example, the first two integrators have dacs with `[1,2]`
   timing, the third has a pair of dacs, one with `[0.5 1]` timing and the
   other with `[1 1.5]` timing, and there is no direct feedback DAC to the
   quantizer
 - `ordering`: A vector specifying which NTF zero-pair to use in each
   resonator. Default is for the zero-pairs to be used in the order specified
   in the NTF.
 - `bp`: A vector specifying which resonator sections are bandpass.
   The default (zeros(...)) is for all sections to be lowpass.
 - `ABCDc`: The loop filter structure, in state-space form.
   If this argument is omitted, ABCDc is constructed according to `form`.
"""
function realizeNTF_ct(ntf, form=:FB, tdac=[0 1], ordering=1:0, bp=Int[], ABCDc=Float64[])
	(ntf_z, ntf_p) = _zpkdata(ntf, 1)
	ntf_z = deepcopy(ntf_z); ntf_p = deepcopy(ntf_p) #Don't modify original
	ntf_z = realfirst(ntf_z); ntf_p = realfirst(ntf_p) #VERIFYME: Expected by realizeNTF... here too?
	ntf_z_orig = deepcopy(ntf_z) #VERIFYME: Doesn't get modified to match legacy code. Intentional or omission??
	order = length(ntf_p)
	(order2, isodd, Δodd) = orderinfo(order)
	if isempty(ordering); ordering = 1:order2; end
	if isempty(bp); bp = zeros(Int, order2); end

	#compensate for limited accuracy of zero calculation
	ntf_z[findall( abs.(ntf_z .- 1) .< eps(1.0)^(1/(1+order)) )] .= 1.0
#ntf_z_orig=ntf_z #???

	multi_tmg = (eltype(tdac)<:Array)
	if multi_tmg
		if size(tdac) != (order+1,)
			error("For multi-array tdac, tdac be a vector of length `order+1`")
		end
		if form != :FB
			error("Currently only supporting form=:FB for multi-array tdac")
		end
	else
		if size(tdac) != (1, 2)
			#TODO: Update message
			error("For non cell array tdac, size(tdac) must be (1, 2)")
		end
	end

	n_direct = 0
	n_extra = 0
	tdac2 = Float64[-1 -1]
	if !multi_tmg
		#Need direct terms for every interval of memory in the DAC
		n_direct = ceil(Int, tdac[2]) - 1
		if ceil(Int, tdac[2]) - floor(Int, tdac[1]) > 1
			n_extra = n_direct-1 #tdac pulse spans a sample point
		else
			n_extra = n_direct
		end
		tdac2 = Float64[ -1 -1; 
			tdac;
			0.5*ones(n_extra,1)*[-1 1] .+ cumsum(ones(n_extra,2), dims=1)
		]
	end

	if isempty(ABCDc)
		ABCDc = zeros(order+1, order+2);
		#Stuff the A portion
		if isodd
			ABCDc[1,1] = real( log( ntf_z[1] ) )
			ABCDc[2,1] = 1
		end
		for i = 1:order2
			n = bp[i]
			i1 = 2*i + Δodd - 1
			zi = 2*ordering[i] + Δodd - 1
			w = abs( angle( ntf_z[zi] ) )
			ABCDc[i1 .+ [0 1 2], i1 .+ [0 1]] = [
				0  -w^2;
				1   0;
				n  1-n;
			]
		end
		ABCDc[1,order+1] = 1
		ABCDc[1,order+2] = -1 #2006.10.02 Changed to -1 to make FF STF have +ve gain at DC
	end

	Ac = ABCDc[1:order, 1:order]
	local Bc, Cc, Dc, tp

	if :FB == form
		Cc = ABCDc[[order+1],1:order]
		if !multi_tmg
			Bc = [eye(order) zeros(order,1)]
			Dc = [zeros(1,order) 1]
			tp = repeat(tdac,order+1,1)
		else	#Assemble tdac2, Bc and Dc
			tdac2 = Float64[-1 -1]
			Bc = ones(order, 0) #0-column
			Bci = [eye(order) zeros(order,1)]
			Dc = ones(1, 0)
			Dci = [zeros(1,order) 1]
			for i in 1:length(tdac)
				tdi = tdac[i]
				if (eltype(tdi)<:Array)
					for j in 1:length(tdi)
						tdj = tdi[j]
						tdac2 = [tdac2; tdj]
						Bc = [Bc Bci[:,i]]
						Dc = [Dc Dci[:,i]]
					end
				elseif !isempty(tdi)
					tdac2 = [tdac2; tdi]
					Bc = [Bc Bci[:,i]]
					Dc = [Dc Dci[:,i]]
				end
			end
			tp = tdac2[2:end,:]
		end
	elseif :FF == form
		Cc = [eye(order); zeros(1,order)]
		Bc = [-1; zeros(order-1,1)]
		Dc = [zeros(order,1); 1]
		tp = tdac #2008-03-24 fix from Ayman Shabra
	else
		error("Sorry, no code for form ", form)
	end

	#Sample the L1 impulse response
	n_imp = ceil(Int, 2*order + maximum(tdac2[:,2]) + 1)
	y = impL1(ntf,n_imp)

	sys_c = _ss( Ac, Bc, Cc, Dc )
	yy = pulse(sys_c, tp, 1.0, 1.0*n_imp, true)
	yy = squeeze(yy)
	#Endow yy with n_extra extra impulses.
	#These will need to be implemented with n_extra extra DACs.
	#!! Note: if t1=int, matlab says pulse(sys) @t1 ~=0
	#!! This code corrects this problem.
	if n_extra>0
		y_right = padb([zeros(1,n_extra+1); eye(n_extra+1)], n_imp+1)
		#Replace the last column in yy with an ordered set of impulses
		yy = [yy[:,1:end-1] y_right[:,end:-1:1]]
	end

	#Solve for the coefficients:
	x = yy\y
	if norm(yy*x - y) > 1e-4
		@warn("Pulse response fit is poor.")
	end
	local Bc2, Dc2
	if :FB == form
		if !multi_tmg
			Bc2 = [ x[1:order] zeros(order,n_extra) ]
			Dc2 = x[order+1:end]'
		else
			BcDc = [Bc; Dc]
			i = findall(BcDc .!= 0)
			BcDc[i] = x
			Bc2 = BcDc[1:end-1,:]
			Dc2 = BcDc[[end],:]
		end
	elseif :FF == form
		Bc2 = [Bc zeros(order,n_extra)]
		Cc = x[1:order]'
		Dc2 = x[order+1:end]'
	else
		error("No code for form \"$form\"")
	end

	Dc1 = 0
	Dc = [Dc1 Dc2]
	Bc1 = [1; zeros(order-1,1)]
	Bc = [Bc1 Bc2]

	#Scale Bc1 for unity STF magnitude at f0
	fz = angle.(ntf_z_orig)/(2π)
	f1 = fz[1]
	ibz = findall(abs.(fz .- f1) .<= abs.(fz .+ f1))
	fz = fz[ibz]
	f0 = mean(fz)
	if minimum(abs.(fz)) < 3*minimum(abs.(fz .- f0))
		f0 = 0
	end
	L0c = _zpk(_ss(Ac,Bc1,Cc,Dc1))
	G0 = evalTFP(L0c, ntf, f0)
	Bc[:,1] = Bc[:,1]/abs.(G0)

	ABCDc = [Ac Bc; Cc Dc]

	return (ABCDc, tdac2)
end

#Last line
