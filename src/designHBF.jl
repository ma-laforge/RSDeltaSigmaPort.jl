#RSDeltaSigma: designHBF() algorithm
#-------------------------------------------------------------------------------
#=IMPORTANT
designHBF6() (a version that doesn't need firpm()) no longer exists in the
source code used in this port. It actually just copies designHBF7().
=#


#==Helper functions
===============================================================================#
"""`F1 = evalF1(f1, z, phi)`

Calculate the values of the F1 filter (tranformed prototype filter) of a
Saramaki HBF at the given points.
"""
function evalF1(f1, z)
	F1 = 0.5
	for i = 1:length(f1)
		F1 = F1 .+ f1[i] * z .^ (2*i-1)
	end
	return F1
end
evalF1(f1, z, phi::Float64) = evalF1(f1, z/phi)

"""`F0 = evalF0(f1, z, phi)`

Calculate the values of the F0 (prototype) filter of a Saramaki HBF
at the given points.
"""
evalF0(f1, z, phi::Float64) =
	evalF1( f1, 0.5*(z + 1 ./ z), phi )


#==frespHBF()
===============================================================================#
function _frespHBF(f, f1, f2, fp, calcripples::Bool)
	Npts = length(f)
	w = 2*pi*f
	z = exp.(j*w)
	cos_w = real.(z)

	n2 = length(f2)
	F2 = zeros(size(w))
	for i = 1:n2
		F2 = F2 + f2[i]*cos.(w*(2*i-1))
	end
	F2 = F2*2
	mag = evalF1(f1, F2)

	pbr = sbr = nothing
	if calcripples
		passband = 1:floor(Int, 2*fp*(Npts-1)+1)
		stopband = (Npts + 1) .- passband
		pbr = maximum( abs.(abs.(mag[passband]) .- 1) )
		sbr = maximum( abs.(mag[stopband]) )
	end
	return (z, mag, F2, pbr, sbr)
end

"""`(mag, pbr, sbr) = frespHBF(f,f1,f2,phi=1.0,fp=0.2,msg="")`

Compute the frequency response, the passband ripple and the stopband ripple 
of a Saramaki HBF. If msg is non-null, a plot is made.
 - `fp` is the passband edge.
 - `phi` is used by designHBF.
"""
function frespHBF(f,f1,f2; phi=1.0, fp=0.2, msg="", calcripples=true)
	if isa(f, Number) && isnan(f)
		f = range(0, stop=.5, length=1024)
	end
	if isa(f1, Vector{CoeffCSD})
		f1 = values(f1)
	end
	if isa(f2, Vector{CoeffCSD})
		f2 = values(f2)
	end
	calcripples = calcripples || isempty(msg)
	(z, mag, F2, pbr, sbr) = _frespHBF(f, f1, f2, fp, calcripples)

	if !isempty(msg)
		F1 = evalF0(f1, z, phi)
		plot = cons(:plot, linlin, nstrips=2, title="Frequency Response", legend=true,
			xaxis=set(label="Normalized Frequency", min=0, max=0.5),
			ystrip1=set(axislabel="Mag", min=0, max=1.1),
			ystrip2=set(axislabel="Mag [dB]", min=-150, max=10),
		)
		plot.title = msg
		w1 = waveform(f, abs.(F1))
		w2 = waveform(f, phi*abs.(F2))
		w3 = waveform(f, abs.(mag))
		w4 = waveform(f, dbv.(mag))
		pbr_str = @sprintf("pbr = %.1e", pbr) #orig. coords: (0.0,-10)
		sbr_str = @sprintf("sbr = %.0fdB", dbv(sbr)) #orig. coords: (0.0,dbv(sbr))

		push!(plot,
			cons(:wfrm, w1, line=set(style=:dash, color=:blue, width=2), label="F1"),
			cons(:wfrm, w2, line=set(style=:dot, color=:red, width=2), label="F2"),
			cons(:wfrm, w3, line=set(style=:solid, color=:green, width=2), label="HBF"),
			cons(:wfrm, w4, line=set(style=:solid, color=:green, width=2), label="HBF", strip=2),
			cons(:atext, pbr_str, y=-10, reloffset=set(x=0.05), align=:tl, strip=2),
			cons(:atext, sbr_str, y=dbv(sbr), reloffset=set(x=.95), align=:br, strip=2),
		)
		displaygui(plot)
	end

	return (mag,pbr,sbr)
end


#==designF1()
===============================================================================#
"""`(f1, zetap, phi) = designF1(delta, fp1)`

Design the F1 sub-filter of a Saramaki halfband filter.
This function is called by designHBF()

 - `f1`: A structure array containing the F1 filter coefficents and
   Their CSD representation.
 - `phi`: The scaling factor for the F2 filter (imbedded in the f1 coeffs.)
"""
function designF1(delta, fp1)
	passband = exp.(4π*j*range(0, stop=fp1, length=10))
	ok = false
	n1_found = 7
	local h
	for n1 in 1:2:7 #Odd values only
		if n1 == 1
			h = [0.5, 0.5]
		else
			h = remez(2*n1, [(0,4*fp1)=>1, (1,1)=>0], Hz=2)
			if !(abs(sum(h)-1) < 1e-3) #remez bug! Use firls instead
				throw("Remez bug! workaround not yet implemented")
				#h = firls(2*n1-1,[0 4*fp1 1-1e-6 1],[1 1 0 0]) #workaround
			end
		end
		fresp = abs.( polyval(h,passband) )
		if maximum( abs.(fresp .- 1) ) <= delta
			n1_found = n1
			ok = true
			break
		end
	end
	n1 = n1_found
	if !ok
		zetap = 1 #Use this as an indication that the function failed.
		return ([], zetap, 0)
	end

	#Transform h(n) to a chebyshev polynomial f1(n)
	#Sum(f1(i)*cos(w)^n)|i=1:n1 + Sum(h(n1+i))*cos(n*w))|i=1:n1, n = 2*i-1;
	w = π*rand(1,n1)
	cos_w = cos.(w)
	A = zeros(n1,length(w))
	B = zeros(1,n1)
	for i in 1:n1
		n = 2*i-1
		A[i,:] = cos_w .^ n
		B = B .+ h[n1+i] * cos.(n*w)
	end
	f1 = B/A
	f1 = f1[:] #Vectorize

	phivecb = []

	#Optimize the quantized version of f1 to maximize the stopband width 
	#( = acos(zetap) )
	zetap = 1
	phi = 0
	f1_saved = nothing
	testPoints = vcat(0, 10 .^ range(-2, stop=0, length=128)) .- 1
	for nsd in 3:8
		f1a = f1; f1b = f1 #First try the unperturbed filter.
		for phia in 1 ./ vcat(1, f1)
			phia = phia / 2^nextpow(2, abs(phia)) #keep phi in (0.5,1]
			#Try a bunch of coefficients in the current neighborhood,
			#shrinking the neighborhood once 10 successive trial values show no
			#improvement.  If 2 successive shrinkages do no good, try a higher nsd.
			count = 0
			nohelp = 0
			neighborhood = .05
			while neighborhood > 1e-5
				phivec = phia .^ collect(1:2:2*n1-1)
				#Matlab Ver. 5 change:
				if isempty(phivecb); phivecb = phivec; end
				f1q = bquantize( f1a .* phivec, nsd )
				f1qval = values(f1q)
				F1 = evalF1( f1qval, testPoints, phia )
				fi = findall( abs.(F1) .> delta )
				zeta = -testPoints[ max(fi[1]-1, 1) ]
				#fprintf(2,'nsd=%d, nbhd= %f, count=%d, zeta = %f, phia=%f\n', ...
				#  nsd, neighborhood, count, zeta, phia )
				if zeta < zetap
					count = 0
					nohelp = 0
					zetap = zeta
					f1b = f1qval
					f1_saved = f1q
					phi = phia
					phivecb = phivec
				else
					count = count + 1
				end
				if count > 10
					count = 0
					neighborhood = neighborhood/2
					nohelp = nohelp + 1
					if nohelp > 2
						break
					end
				end
				f1a = f1b ./ phivecb + neighborhood*(rand(size(f1b)...) .- 0.5)
				phia = phia .+ neighborhood*(rand()-0.5)
			end #while neighborhood ...
			if zetap < 1 #Found a filter with adequate attn.
				break
			end
		end #for phia ...
		if zetap < 1 #Found a filter with adequate attn.
			break
		end
	end #for nsd ...

	return (f1_saved, zetap, phi)
end


#==designF2()
===============================================================================#
"""`f2 = designF2(fp,zetap,phi)`

Design the F2 sub-filter of a Saramaki halfband filter.
This function is called by designHBF().

# subfilter design:
   `1 - delta2' < |F2/phi| < 1`  for `f` in `[0 fp]`
   `-1 < |F2/phi| < -1 + delta2` for `f` in `[0.5-fp, 0.5]`
   `1-delta2 = (1-delta2)/(1+delta2)`
"""
function designF2(fp, zetap, phi)
	delta2 = (1-zetap)/(1+zetap)
	#delta2p = 1 - (1-delta2)/(1+delta2)

	#determine the minimum order required by the filter
	passband = exp.(j*range(0, stop=4π*fp, length=10))
	nsub_found = 17
	for nsub in 3:2:17
		h2 = remez(nsub+1, [(0,4*fp)=>1, (1,1)=>0], Hz=2)
		mag = abs.( polyval(h2, passband) )
		if maximum(abs.(mag .- 1)) < delta2
			nsub_found = nsub
			break
		end
	end
	n2min = (nsub_found+1)/2

	#Search all n2,nsd pairs, in order of the product n2*(nsd+1)
	#allowing fp to be a variable?
	success = false
	nsdmin = 3; nsdmax = 6
	for product in (nsdmin+1)*n2min:(nsdmax+1)*n2min
		for nsd in nsdmin:nsdmax
			n2 = product/(nsd+1)
			if floor(n2) != n2 #Only take integer n2,nsd pairs
				break
			end
			n2 = floor(Int, n2)
			nsub = 2*n2-1
			#Could try a bunch of fp values
			#fprintf(2,'designF2: Trying (n2,nsd2,fp)=(%2d,%2d,%6.4f)\n',n2,nsd,fp);
			h2 = remez(nsub+1, [(0,4*fp)=>1, (1,1)=>0], Hz=2)
			h2 =  h2/(phi*(1+delta2)) #Adjust the coefficients.
			f2 = bquantize( h2[n2+1:nsub+1], nsd )
			f2val = values(f2)
			h2 = (1+delta2)*phi*vcat(f2val[n2:-1:1], f2val)
			mag = abs.( polyval(h2, passband) )
			if maximum(abs.(mag .- 1)) < delta2
				success = true
				break
			end
		end
		if success
			break
		end
	end

	if !success
		f2 = []
		q2 = []
	end
	return f2
end


#==designHBF()
===============================================================================#
"""`(f1,f2,info)=designHBF(fp(=0.2),delta=1e-5,debug=false)`

Design a half-band filter which can be realized without general multipliers.
The filter is a composition of a prototype and sub- filter.

# Input
 - `fp`: The normalized cutoff frequency of the filter. Due to the
   symmetry imposed by a HBF, the stopband begins at `0.5-fp`.
 - `delta`: The absolute value of the deviation of the frequency response from 
   the ideal values of 1 in the passband and 0 in the stopband.

# Output
 - `f1`,`f2`: The coefficients of the prototype and sub-filters
   and their canonical-signed digit (csd) representation.
 - `info`: A vector containing the following data (only set when debug=1):
 - -->`complexity`: The number of additions per output sample.
 - -->`n1,n2`: The length of the f1 and f2 vectors.
 - -->`sbr`: The achieved stob-band attenuation (dB).
 - -->`phi`: The scaling factor for the F2 filter.
"""
function designHBF(fp::Float64=0.2; delta::Float64=1e-5, debug::Bool=false)
	f1_saved = f2_saved = nothing #Declare variables
	phi_saved = 0.0
#=To Do:
	Clean up the code a bit more, esp. wrt the use of the struct. arrays.
	Use the phi variable to cut down on the number of adders in F2.
	Apply a simulated annealing/genetic optimization alg instead
	of the ad hoc one I have now.
=#

	#Try several different values for the fp1 parameter.
	#The best values are usually around .04
	#Surrender if 3 successive attempts yield progressively greater complexity.
	lowest_complexity = Inf; prev_complexity = Inf;
	local worse
	for fp1 in [.03, .035, .025, .040, .020, .045, .015, .05]
		failed = false
		(f1, zetap, phi) = designF1(delta, fp1)
		if zetap == 1 #designF1 failed
			failed = true
			if debug
				@warn("designF1 failed at fp1=$fp1\n")
			end
		end
		if !failed
			f2 = designF2(fp, zetap, phi)
			n1 = length(f1);	n2 = length(f2)
			if n2 == 0 #designF2 failed
				failed = true
				if debug
					@warn("designF2 failed when zetap=$zetap, phi=$phi\n")
				end
			end
		end
		if !failed
			#complexity(+ performance) = the number of two-input adders (+ sbr)
			complexity = sum(csdsize2(f1)) + (2*n1-1)*(n2+sum(csdsize2(f2))-1) #VERIFYME
			msg = ""

			if debug
				msg = @sprintf("%d adders: n1=%d, n2=%d, (fp1=%.2f, zetap=%.3f, phi=%4.2f)",
					complexity, n1, n2, fp1, zetap, phi
				)
			end
			(fresp, pbr, sbr) = frespHBF(NaN, f1, f2, phi=phi, fp=fp, msg=msg)
			if pbr <= delta && sbr <= delta
				complexity = complexity + sbr
				if complexity < prev_complexity
					worse = 0
					if complexity < lowest_complexity 
						lowest_complexity = complexity
						f1_saved = f1;	f2_saved = f2
						phi_saved = phi
						if debug; println(msg); end
					end
				else
					worse = worse + 1
					if worse > 2; break; end
				end
				prev_complexity = complexity
			end #if pbr <= delta
		end #!failed
	end #for fp1

	info=nothing
	if isinf(lowest_complexity)
		@warn("designHBF: Unable to meet the design requirements.")
	else
		complexity = floor(lowest_complexity)
		msg = "Final Design: $complexity adders"
		(junk, pbr, sbr) = frespHBF(NaN, f1_saved, f2_saved, phi=phi_saved, fp=fp, msg=msg)
		n1 = length(f1_saved); n2 = length(f2_saved)
		if debug
			msg *= @sprintf(" (%d,%d,%.0fdB)", n1, n2, dbv(sbr))
			println(msg)
		end
		info = (complexity=complexity, n1=n1, n2=n2, sbr=dbv(sbr), phi=phi_saved)
	end

	return (f1_saved, f2_saved, info)
end

#Last line
