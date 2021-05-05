# Demonstrate the simulateMS function
using RSDeltaSigmaPort
using RSDeltaSigmaPort.EasyPlot #set, cons
import RSDeltaSigmaPort: fft
import RSDeltaSigmaPort: BoundingBox
import Printf: @sprintf
j=im


#==Demo5Module: Module to define types, etc
===============================================================================#
module Demo5Module

struct SimConfig
	A::Float64
	f::Float64
	mtf
	dither::Float64
#	dw
	id::String
end

struct SimResult
	fin::Int #Input frequency
	sv::Array
	Svv::Array
	Sdd::Array
	leg::String
end

end #Demo5Module

import .Demo5Module: SimConfig, SimResult


#==
===============================================================================#
show_usage = false
OSR = 25
M = 16
N = 2^14
sigma_d = 0.01 #1% mismatch

dsm2b = RealDSM(order=2, OSR=OSR, M=M, opt=1, Hinf=2) #Second-order shaping
dsm4 = RealDSM(order=4, OSR=round(Int, OSR*0.9), M=M, opt=1, Hinf=1.3) #Fourth-order shaping
dsm6 = RealDSM(order=6, OSR=OSR, M=M, opt=1, Hinf=4)

mtf1 = _zpk(1,0,1,1)         #First-order shaping
mtf2 = _zpk([ 1 1 ], [ 0.3 0.3 ], 1, 1) #Second-order shaping
mtf2b = synthesizeNTF(dsm2b) #Second-order shaping
mtf4 = synthesizeNTF(dsm4)   #Fourth-order shaping

cases = [
	SimConfig(undbv(-3),  0.01,    [],   0, "Thermometer"),
	SimConfig(undbv(-3),  0.01,  mtf1,   0, "Rotation"),
	SimConfig(undbv(-3),  0.01,  mtf2,   0, "2^{nd}-order"),
	SimConfig(undbv(-30), 0.01,  mtf1,   0, "Rotation"),
	SimConfig(undbv(-30), 0.01,  mtf1, 0.5, "Rot + dither"),
	SimConfig(undbv(-3),  0.01, mtf2b,   0, "2^{nd}-order with zero"),
	SimConfig(undbv(-3),  0.01,  mtf4,   0, "4^{th}-order"),
]
comparisons = [ [1,2], [2,3], [4,5], [3,6], [6,7] ]

if false #Debug
	cases = cases[1:2]
	comparisons = comparisons[1:1]
end

nconfig = length(cases)
println("\nNumber of configurations to simulate: $nconfig"); flush(stdout); flush(stderr)
#-------------------------------------------------------------------------------

NTF = synthesizeNTF(dsm6)
window = ds_hann(N)'/(M*N/8)
windowM = repeat(window, M)

println()
@info("Running simulations..."); flush(stdout); flush(stderr)
#-------------------------------------------------------------------------------
simresults = SimResult[] #Simulation results
for cfg in cases
	local u #Avoid complaints wrt "global u"
	local sv
	println("\t", cfg.id); flush(stdout); flush(stderr)
	fin = round(Int, cfg.f*N)
	inband = setdiff(2:1 .+ ceil(Int, 0.5*N/OSR), 1 .+ [0 1 fin .+ [-1 0 1]])
	w = (2π/N)*fin
	u = M*cfg.A*sin.(w*(0:N-1))
	v = simulateDSM(u, NTF, nlev=M+1).v #M unit elements requires an (M+1)-level quant.
	data = applywnd(v, window)
	Svv = abs.(fft(applywnd(v, window))) .^ 2
	if isa(cfg.mtf, Array) && isempty(cfg.mtf)
		sv = ds_therm(v, M)
	else
		sv = simulateMS(v, M=M, mtf=cfg.mtf, d=cfg.dither).sv
	end
	Sdd = sigma_d^2 * sum( abs.(fft(applywnd(sv, windowM),2)) .^ 2 , dims=1)
	mnp = sum(Sdd[inband])/1.5
	leg = @sprintf("%s (MNP= %.0f dBFS)", cfg.id, dbp(mnp))
	push!(simresults, SimResult(fin, sv, Svv, Sdd, leg))
end
println("\tdone.")

println()
@info("Plotting results..."); flush(stdout); flush(stderr)
#-------------------------------------------------------------------------------
for case_nums in comparisons
	local plot #Avoid complaints wrt "global plot"
	nc = length(case_nums)
	#_cases = cases[case_nums]
	_results = simresults[case_nums]
	id1 = _results[1].leg; id2 = _results[2].leg #Only print out first 2
	println("Mismatch-Shaping Unit-Element DAC")
	println("\tComparing $id1 vs. $id2")
	if show_usage
		pcoll = cons(:plot_collection, ncolumns=1, title="Element Usage")
		T = 25
		for resi in _results
			plot = plotUsage(resi.sv[:,1:T])
			push!(pcoll, plot)
		end
		displaygui(pcoll)
	end

	colors = [:blue, :magenta, :red, :green, :orange]
	plot = plotSpectrum()
#		plot.title = "Error Spectra"
	set(plot,
		xyaxes=set(xmin=1e-3, xmax=0.5, ymin=-140, ymax=-50),
		labels=set(yaxis="Error PSD"),
	)
	for (i, resi) in enumerate(_results)
		plotSpectrum!(plot, sqrt.(resi.Sdd), resi.fin, lw=2, color=colors[i], id=resi.leg, n=4)
	end
	resi = _results[end]
	A = dbp(resi.Svv[resi.fin])
	plotSpectrum!(plot, sqrt.(resi.Svv), resi.fin, color=:green, n=5)
	lwfrm = waveform([1e-3, 0.5/OSR], -140*[1, 1])
	nbw_str = @sprintf("NBW = %.1e", 1.5/N) #orig. coords: (0.5, -140)
	push!(plot,
		cons(:wfrm, lwfrm, line=set(style=:solid, color=:black, width=8)),
		cons(:atext, nbw_str, y=-140, reloffset=set(x=0.95), align=:br),
	)
	plot.title = @sprintf("A = %.0fdBFS", A)
	displaygui(plot)
end

if true #Quadrature example
	println()
	@warn("TODO: Quadrature example")
end


:END_OF_DEMO
