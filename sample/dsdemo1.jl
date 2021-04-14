# Demonstrate NTF synthesis (synthesizeNTF)
using RSDeltaSigmaPort
j=im


#==Baseband modulator
===============================================================================#
println("\n*** 5th order, 2-level, baseband modulator")
OSR = 32

@info("Synthesizing NTF without zero-optimization..."); flush(stdout); flush(stderr)
#-------------------------------------------------------------------------------
NTF_noopt = synthesizeNTF(5, OSR, opt=0)
println("\tdone.")

plot = plotNTF(NTF_noopt, OSR, color=:blue)
plot.title = "5th-Order Modulator (No Zero-Optimization)"
saveimage(:png, "dsdemo1_o5_noopt.png", plot, AR=2/1, width=900)
displaygui(plot)

@info("Synthesizing NTF with optimized zeros..."); flush(stdout); flush(stderr)
#-------------------------------------------------------------------------------
NTF_opt = synthesizeNTF(5, OSR, opt=1)
println("\tdone.")

plot = plotNTF(NTF_opt, OSR, color=:red)
plot.title = "5th-Order Modulator (Optimized Zeros)"
saveimage(:png, "dsdemo1_o5_zopt.png", plot, AR=2/1, width=900)
displaygui(plot)

@info("Plotting NTF comparison (overlay results)")
#-------------------------------------------------------------------------------
plot = plotNTF(NTF_noopt, OSR, color=:blue)
plot = plotNTF!(plot, NTF_opt, OSR, color=:red)
plot.title = "5th-Order Modulator (Optimized Zeros - Overlay)"
saveimage(:png, "dsdemo1_o5_cmp.png", plot, AR=2/1, width=900)
displaygui(plot)


#==Bandpass modulator
===============================================================================#
println("\n*** 8th order, 2-level, bandpass modulator")
OSR = 64
order = 8
f0 = 0.125 #fs/8

function calcSTF(order, OSR, NTF, f0)
	G = _zpk(zeros(array_round(order/2)),NTF.p,1,1)
	G.k = 1/abs(evalTF(G,exp(2π*j*f0)))
	return G
end

@info("Synthesizing NTF..."); flush(stdout); flush(stderr)
#-------------------------------------------------------------------------------
NTF = synthesizeNTF(order, OSR, opt=2, f0=f0)
println("\tdone.")

@info("Plotting NTF")
#-------------------------------------------------------------------------------
STF = calcSTF(order, OSR, NTF, f0)
plot = plotNTF(NTF, OSR, color=:blue, f0=f0, STF=STF)
plot.title = "8th-Order Bandpass Modulator"
saveimage(:png, "dsdemo1_o8_bp.png", plot, AR=2/1, width=900)
displaygui(plot)

:END_OF_DEMO
