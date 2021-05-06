#Demonstrate Saramaki half-band filter design
using RSDeltaSigmaPort
using RSDeltaSigmaPort.EasyPlot #set, cons
using RSDeltaSigmaPort: csdsize2, linlin
import RSDeltaSigmaPort: BoundingBox
import Printf: @sprintf
j=im

use_canned_example = true
println("\n*** Half-band filter design (ALPHA VERSION)")
f1=f2=0 #Define variables
title=nothing
if use_canned_example
	(f1, f2) = exampleHBF(2)
else
	if !RSDeltaSigmaPort.FIRPM_AVAIL
		throw("firpm() not available. Cannot run designHBF().")
	end
	fp = 0.9*0.25
	delta = undbv( -100 )
	title = "designHBF Iterations"
	(f1, f2, info) = designHBF(fp,delta,1)
end
n1 = length(f1); n2 = length(f2)
complexity = sum(csdsize2(f1)) + (2*n1-1)*(n2+sum(csdsize2(f2))-1) #VERIFYME
@show complexity

println()
@info("Running simulations..."); flush(stdout); flush(stderr)
#-------------------------------------------------------------------------------
#interleave the even and odd decimated impulse responses
Nimp = 2^11
imp = simulateHBF(vcat(1, zeros(Nimp-1)),f1,f2)
mag = abs.(fft(imp))
mag = mag[1:array_round(end/2+1)]
t = range(0, stop=0.5, length=length(mag))
println("\tdone.")

println()
@info("Plotting resulting filter"); flush(stdout); flush(stderr)
#-------------------------------------------------------------------------------
plot = cons(:plot, linlin, title="designHBF() result", legend=false,
	xyaxes=set(xmin=0, xmax=0.5, ymin=-150, ymax=10),
	labels=set(xaxis="Normalized Frequency", yaxis="|H| [dB]"),
)
magdBw = waveform(t, dbv.(mag))
push!(plot,
	cons(:wfrm, magdBw, line=set(style=:solid, color=:blue, width=2), label=""),
)
displaygui(plot)

:END_OF_DEMO
