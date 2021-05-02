#RSDeltaSigma: mapCtoD (State space utility)
#-------------------------------------------------------------------------------


#==Helper functions
===============================================================================#
function B2formula(Ac, t1, t2, B2)
	if t1==0 && t2==0
		term = B2
		return term
	end
	n = size(Ac, 1)
	tmp = eye(n) - exp(-Ac)
	if cond(tmp) < 1e6
		term = ( (exp(-Ac*t1) - exp(-Ac*t2)) * inv(tmp)) * B2
		return term
	end
	#Numerical trouble. Perturb slightly and check the result
	ntry = 0
	k = sqrt(eps)
	Ac0 = Ac
	#Make this operation repeatable
	#rng_state = rng; rng('default') #Original code
	Random.seed!(0)
	while ntry <2
		Ac = Ac0 + k*rand(n,n)
		tmp = eye(n) - exp(-Ac)
		if cond(tmp) < 1/sqrt(eps)
			ntry = ntry+1
			if ntry == 1
				term = ( (exp(-Ac*t1) - exp(-Ac*t2)) * inv(tmp)) * B2
			else
				term1 = ( (exp(-Ac*t1) - exp(-Ac*t2)) * inv(tmp)) * B2
			end
		end
		k = k*sqrt(2)
	end
	Random.seed!() #original code restores state: rng(rng_state)
	if norm(term1-term) > 1e-3
		@warn("Inaccurate calculation in mapCtoD.")
	end
	return term
end


"""`(sys, Gp) = mapCtoD(sys_c, t=[0 1], f0=0, calcGp=false)`

Map a MIMO continuous-time system to a SIMO discrete-time equivalent.
The criterion for equivalence is that the sampled pulse response
of the CT system must be identical to the impulse response of the DT system.
i.e. If yc is the output of the CT system with an input vc taken
from a set of DACs fed with a single DT input v, then y, the output
of the equivalent DT system with input v satisfies
y(n) = yc(n-) for integer n. The DACs are characterized by
rectangular impulse responses with edge times specified in the t matrix.

# Input
 - `sys_c`: The LTI description of the CT system.
 - `t`: The edge times of the DAC pulse used to make CT waveforms
   from DT inputs. Each row corresponds to one of the system
   inputs; `[-1 -1]` denotes a CT input. The default is `[0 1]`,
   for all inputs except the first.
 - `f0`: The (normalized) frequency at which the Gp filters' gains are
   to be set to unity. Default 0 (DC).

# Output
 - `sys`: The LTI description for the DT equivalent
 - `Gp`: The mixed CT/DT prefilters which form the samples
   fed to each state for the CT inputs.

# Reference
 - R. Schreier and B. Zhang, "Delta-sigma modulators employing continuous-time
   circuitry," IEEE Transactions on Circuits and Systems I, vol. 43, no. 4,
   pp. 324-332, April 1996.
"""
function mapCtoD(sys_c; t=[0 1], f0::Float64=0.0, calcGp::Bool = false)
	ni = size(sys_c.B, 2)

	if size(t) == (1,2) && (ni > 1) #VERIFYME
		t = [-1 -1; ones(ni-1,1)*t]
	end
	if size(t) != (ni, 2) #VERIFYME
		error("The t argument has the wrong dimensions.")
	end

	di = ones(1,ni)
	for i in 1:ni
		if t[i,:]==[-1 -1]
			di[i] = 0
		end
	end
	di = (di .!= 0) #Convert to logical array.
	i!di = findall(.!di); idi = findall(di)

	#c2d assumes t1=0, t2=1.
	#Also c2d often complains about poor scaling and can even produce
	#incorrect results.
	NEEDS_DELAY = false #VERSIONSWITCH
	if NEEDS_DELAY
		#set(sys_c,'InputDelay',0)
	end
	sys = _c2d(sys_c, 1, :zoh)
	Ac = sys_c.A
	Bc1 = sys_c.B[:, i!di]
	A = sys.A; B = sys.B; C = sys.C; D = sys.D
	B1 = B[:, i!di]; D1 = D[:, i!di]

	#Examine the discrete-time inputs to see how big the
	#augmented matrices need to be.
	n = size(A,1)
	t2 = ceil(t[idi,2])
	esn = (t2 .== t[idi,2]) .& (D(1,idi) .!= 0)' #extra states needed?
	np = n + maximum(t2 .- 1 + esn)

	#Augment A to np x np, B to np x 1, C to 1 x np.
	Ap = padb( padr(A,np), np )
	for i in n+2:np
		Ap[i, i-1] = 1
	end
	Bp = zeros(np, 1)
	if np>n
		Bp[n+1] = 1
	end
	Cp = padr(C, np)
	Dp = 0

	#Add in the contributions from each DAC
	for i in findall(di)
		t1 = t[i,1]
		t2 = t[i,2]
		B2 = B[:,i]
		D2 = D[:,i]
		if t1==0 && t2==1 && all(D2 .== 0) #No fancy stuff necessary
			Bp = Bp + padb(B2,np)
		else
			n1 = floor(Int, t1); n2 = ceil(Int, t2)-n1-1
			t1 = t1-n1;	t2 = t2-n2-n1
			if t2==1 && all(D2 .!= 0)
				n2 = n2+1
				extraStateNeeded = true
			else
				extraStateNeeded = false
			end
			nt = n+n1+n2
			if n2>0 && t2!=1
				Ap[1:n,nt] = Ap[1:n,nt] + B2formula(Ac, 0, t2, B2)
			end
			local Btmp
			if n2>0 #pulse extends to the next period
				Btmp = B2formula(Ac,t1,1,B2)
			else    #pulse ends in this period
				Btmp = B2formula(Ac,t1,t2,B2)
			end
			if n1>0
				Ap[1:n,n+n1] = Ap[1:n,n+n1] + Btmp
			else
				Bp = Bp + padb(Btmp,np)
			end
			if n2>0
				Cp = Cp + padr([zeros(size(D2,1),n+n1) D2*ones(1,n2)], np)
			end
		end
	end
	sys = _ss(Ap, Bp, Cp, Dp, 1)

	Gp = nothing
	if any(di .== 0)
		if calcGp
			#Compute the prefilters and add in the CT feed-ins.
			#Gp = inv(sI-Ac) * (zI-A)/z *Bc1
			(n, m) = size(Bc1)
			Gp = fill(CTDTPrefilters[], n,m)
			ztf = Array{Any}(nothing, size(Bc1)) #!!Make this like stf: an array of zpk objects
			#Compute the z-domain portions of the filters
			ABc1 = A*Bc1
			for h in 1:m
				for i in 1:n
					if Bc1[i,h]==0
						ztf[i,h] = _zpk([],0,-ABc1[i,h],1)
					else
						ztf[i,h] = _zpk(ABc1[i,h]/Bc1[i,h],0,Bc1[i,h],1)
					end
				end
			end
			#Compute the s-domain portions of each of the filters
			stf = _zpk(_ss(Ac,eye(n),eye(n),zeros(n,n))) #Doesn't do pole-zero cancelation
			for h in 1:m
				for i in 1:n
					for j=1:n
						if stf.k[i,j]!=0 && ztf[j].k!=0 #A non-zero term
							push!(Gp[i,h], CTDTPrefilters(stf(i,j), ztf[j,h]))
						end
					end
				end
			end
			if f0 != 0 #Need to correct the gain terms calculated by c2d
				#B1 = gains of Gp @f0;
				for h in 1:m
					for i in 1:n
						B1[i,h] = evalMixedTF(Gp[i,h],f0)
						#abs() used because ss() whines if B has complex entries...
						#This is clearly incorrect.
						#I've fudged the complex stuff by including a sign....
						B1[i,h] = abs(B1[i,h]) * sign(real(B1[i,h]))
						if abs(B1[i,h]) < 1e-9
							B1[i,h] = 1e-9 #This prevents NaN in line 174
						end
					end
				end
			end
			#Adjust the gains of the pre-filters
			for h in 1:m
				for i in 1:n
					for j in 1:length(Gp[i,h])
						Gp[i,h][j].Hs.k = Gp[i,h][j].Hs.k/B1[i,h] #This is line 174
					end
				end
			end
			sys.B = [padb(B1,np) sys.B]; sys.D = [D1 sys.D]
		else #Cheat and just dublicate the B1 terms
			sys.B = [padb(Bc1,np) sys.B]; sys.D = [D1 sys.D]
		end
	end

	return (sys, Gp)
end


#Last line
