#RSDeltaSigmaPort: Generate plots for noise transfer functions
#-------------------------------------------------------------------------------


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
"""`axis_handle = documentNTF(ntf|ABCD|mod_struct,osr=64,f0=0,quadrature=0,sizes,plotFreqRespOnly=1)`

The first argument is either the NTF, ABCD matrix or a struct containing 
ntf, osr=64, f0=0, quadrature=0 and optionally stf. 
If the first argument is a struct, then no other arguments should be supplied.
If the first argument is ABCD, the stf is also plotted.
"""
function documentNTF(arg1,osr,f0,quadrature,sizes,plotFreqRespOnly)
	@warn("documentNTF not implemented (see DocumentNTF.m)")
	return #axis_handle
end

#Last line
