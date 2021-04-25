#RSDeltaSigmaPort: Generate plots for noise transfer functions
#-------------------------------------------------------------------------------


#==Vector generators
===============================================================================#
"""`𝑓 = ds_freq(OSR=64, f0=0, quadrature=false)`

Frequency vector suitable for plotting the frequency response of an NTF.
"""
function ds_freq(;OSR::Int=64, f0::Float64=0, quadrature::Bool=false)
	local f_left, f_special
	if quadrature
		f_left = -0.5
		f_special = [f0 -f0]
	else
		f_left = 0
		f_special = f0
	end

	𝑓 = collect(range(f_left, stop=0.5, length=100))
	#Use finer spacing in the vicinity of the passband
	for fx = f_special
		f1 = max(f_left, fx-1/OSR)
		f2 = min(0.5,fx+2/OSR)
		deleteat!(𝑓, (𝑓 .<= f2) .& (𝑓 .>= f1))
		𝑓 = sort( vcat(𝑓, range(f1, stop=f2, length=100)) )
	end
	return 𝑓
end


#==NTF plots
===============================================================================#
function plotNTF(isbandpass)
	dBmin = -100 #Limit dB window to finite value
	𝑓min_log = 10^-2 #Limit min 𝑓 for log(𝑓-norm) plots
	nogrid = set(vmajor = false, vminor = false, hmajor=false, hminor=false)
	vgridonly = set(vmajor = true, vminor = true, hmajor=false, hminor=false)

	plot_PZ = plotPZ()
	plot_PZ.title = "NTF Poles and Zeros"

	plot_NTFlin = cons(:plot, legend=false,
		title = "NTF Magnitude Response",
		xyaxes=set(xscale=:lin, yscale=:lin, ymin=dBmin),
		labels=set(xaxis="Normalized frequency (1→f_s)", yaxis="dB"),
	)

	plot_NTFlog = cons(:plot, legend=false, grid = vgridonly,
#		title = "NTF Magnitude Response",
		xyaxes=set(xscale=:log, yscale=:lin, xmin=𝑓min_log, ymin=dBmin),
		labels=set(xaxis="Normalized frequency (1→f_B)", yaxis="dB"),
	)
	plot_NTFzoom = cons(:plot, legend=false, grid = vgridonly,
#		title = "NTF Magnitude Response",
		xyaxes=set(xscale=:lin, yscale=:lin, xmin=-0.6, ymax=0.6, ymin=dBmin),
		labels=set(xaxis="Normalized frequency offset", yaxis="dB"),
	)
	if !isbandpass; plot_NTFzoom = plot_NTFlog; end

	pcoll = push!(cons(:plot_collection, title="Sample Plot"),
		plot_PZ, plot_NTFlin, plot_NTFzoom
	)

	pcoll.bblist = [
		BoundingBox(0, 0.5, 0, 1), #Pole-zero
		BoundingBox(0.5, 1, 0, 0.5), #NTF linf
		BoundingBox(0.5, 1, 0.5, 1), #NTF zoom
	]

	return pcoll
end

function plotNTF!(pcoll, NTF::ZPKData, OSR::Int; color=:blue, f0::Real=0, STF=nothing)
	dfltline = cons(:a, line=set(style=:solid, color=color, width=2))
	gainline = cons(:a, line=set(style=:solid, color=:red, width=2))
	markerline = cons(:a, line=set(style=:dash, width=2.5))
	plot_PZ, plot_NTFlin, plot_NTFzoom = pcoll.plotlist
	isbandpass = f0!=0

	#Plot poles and zeros
	plotPZ!(plot_PZ, NTF, color=color)

	#Plot full-range NTF
	𝑓 = vcat(range(0, stop=0.75/OSR, length=100), range(0.75/OSR, stop=0.5, length=100))
	if isbandpass
		f1 = f0-1/(2*OSR); f2 = f0+1/(2*OSR)
		𝑓 = vcat(range(0, stop=f1, length=50), range(f1, stop=f2, length=100), range(f2, stop=0.5, length=50))
	end
	z = exp.(2j*π*𝑓)
	magNTF = dbv(evalTF(NTF, z))
	push!(plot_NTFlin, 
		cons(:wfrm, DataF1(𝑓, magNTF), dfltline, label="|NTF|")
	)
	if STF != nothing
		plot_NTFlin.title = "NTF/STF Magnitude Response"
		magSTF = dbv(evalTF(STF,z))
		push!(plot_NTFlin, 
			cons(:wfrm, DataF1(𝑓, magSTF), gainline, label="|STF|")
		)
	end


	#Plot "zoomed-in" NTF
	𝑓start = 0.01
	𝑓 = collect(range(𝑓start, stop=1.2, length=200)/(2*OSR))
	∫rng = [0, 0.5/OSR]
	if isbandpass
		f1 = f0-0.3/OSR; f2 = f0+0.3/OSR
		𝑓 = collect(range(f1, stop=f2, length=100))
		∫rng = f0 .+ [-.25, .25]/OSR
	end
	z = exp.(2j*π*𝑓)
	𝑓norm = (𝑓 .- f0)*(2*OSR)
	magNTF = dbv(evalTF(NTF, z))
	σ_NTF = dbv(rmsGain(NTF, ∫rng[1], ∫rng[2]))
	rmsgain_str = @sprintf("RMS Gain = %5.0fdB", σ_NTF)
	push!(plot_NTFzoom, 
		cons(:wfrm, DataF1(𝑓norm, magNTF), dfltline, label="|NTF|"),
		cons(:atext, rmsgain_str, y=σ_NTF, offset=set(y=3), reloffset=set(x=0.5), align=:bc),
		cons(:hmarker, σ_NTF, markerline, strip=2),
	)

	return pcoll
end

plotNTF(NTF::ZPKData, args...; f0=0, kwargs...) =
	plotNTF!(plotNTF(f0!=0), NTF::ZPKData, args...; f0=f0, kwargs...)


#==documentNTF
===============================================================================#
"""`plot = documentNTF(NTF|ABCD|mod_struct; OSR=64, f0=0, quadrature=false, frespOnly=true, sizes=nothing)`

The first argument is either the NTF, ABCD matrix or a struct containing 
NTF, OSR=64, f0=0, quadrature=0 and optionally stf. 
If the first argument is a struct, then no other arguments should be supplied.
If the first argument is ABCD, the stf is also plotted.
"""
function documentNTF(arg1; OSR::Int=64, f0::Float64=0.0, quadrature::Bool=false, frespOnly::Bool=true, sizes=nothing)
	STF=nothing
	local NTF
	if isa(arg1, ZPKData)
		NTF = arg1
	else #if isa(arg1, Array) #Assume ABCD matrix
		ABCD = arg1
		NTF, STF = calculateTF(ABCD)
	end

	logplot = f0==0
	if quadrature
   	f_left = -0.5
   elseif logplot
   	f_left = 1e-3
   else
   	f_left = 0.0
	end

	pcoll = cons(:plot_collection, ncolumns=2, title="")
	if !frespOnly
		#Plot poles and zeros
		plot = plotPZ(NTF, color=:blue)
		push!(pcoll, plot)

		pcoll.bblist = [ #Format plot locations
			BoundingBox(0, 0.5, 0, 1), #Pole-zero plot
			BoundingBox(0.5, 1, 0, 1), #Frequency spectrum
		]
	end

	#Frequency response
	𝑓 = ds_freq(OSR=OSR, f0=f0, quadrature=quadrature)
	z = exp.(2π*j*𝑓)
	H = dbv.(evalTF(NTF, z))

	axes = logplot ? loglin : linlin
	plot = cons(:plot, axes, title="Frequency Response", legend=false,
		labels=set(xaxis="Normalized Frequency", yaxis=""),
		xyaxes=set(xmin=f_left, xmax=0.5, ymin=-100, ymax=15),
	)
	push!(pcoll, plot)

	H = waveform(𝑓, H)
	push!(plot,
		cons(:wfrm, H, line=set(style=:solid, color=:blue, width=2), label=""),
	)

	if !isnothing(STF)
#		plot.title = "NTF and STF"
		G = dbv.(evalTF(STF,z))
		G = waveform(𝑓, G)
		push!(plot,
			cons(:wfrm, G, line=set(style=:solid, color=:magenta, width=2), label=""),
		)
	end
	f1, f2 = ds_f1f2(OSR=OSR, f0=f0, iscomplex=quadrature)
	NG0 = dbv.(rmsGain(NTF, f1, f2))
	NG0w = waveform([max(f1,f_left), f2], NG0*[1,1])
	push!(plot,
		cons(:wfrm, NG0w, line=set(style=:dash, color=:black, width=2), label=""),
	)

	tstr = @sprintf("  %.0fdB", NG0)
	if f0==0 && logplot
		a = cons(:atext, tstr, x=sqrt(f_left*0.5/OSR), y=NG0+2, align=:bc)
	else
		a = cons(:atext, tstr, x=f2, y=NG0-1, align=:cl)
	end
	push!(plot, a)
	#set(h, 'FontSize',sizes.fs, 'FontWeight',sizes.fw);
	msg = @sprintf(" Inf-norm of H = %.2f\n 2-norm of H = %.2f", infnorm(NTF)[1], rmsGain(NTF,0,1));
	if f0<0.25
		a = cons(:atext, msg, x=0.48, y=0.0, align=:tr)
	else
		a = cons(:atext, msg, x=f_left, y=0.0, align=:tl)
	end
	push!(plot, a)
	#set(h, 'FontSize',sizes.fs, 'FontWeight',sizes.fw);
	if quadrature
		ING0 = dbv(rmsGain(ntf,-f1,-f2))
		ING0w = waveform(-1*[f1, f2], ING0*[1, 1])
		tstr = @sprintf("%.0fdB", ING0)
		push!(plot,
			cons(:wfrm, NG0w, line=set(style=:solid, color=:black, width=2), label=""),
			cons(:atext, tstr, x=-f0, y=ING0+1, align=:bc),
		)
	end
	#set(h, 'FontSize',sizes.fs, 'FontWeight',sizes.fw);

	return pcoll
end

"""`documentNTF(dsm, SYS=nothing; frespOnly=true, sizes=nothing)`
"""
function documentNTF(dsm::AbstractDSM, SYS=nothing; frespOnly::Bool=true, sizes=nothing)
	isnothing(SYS) && (SYS = synthesizeNTF(SYS))
	return documentNTF(SYS; OSR=dsm.OSR, f0=dsm.f0, quadrature=false, frespOnly=frespOnly, sizes=sizes)
end

#Last line
