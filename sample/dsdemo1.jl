# Demonstrate synthesizeNTF
using RSDeltaSigmaPort
j=im

#==5th-order modulator
===============================================================================#
println("\t\tNTF Synthesis-- 5th-Order Modulator\n")
OSR = 32

#Basic
NTF5_noopt = synthesizeNTF(5, OSR, opt=0)
plot = plotNTF(NTF5_noopt, OSR, color=:blue)
plot.title = "5th-Order Modulator"
saveimage(:png, "dsdemo1_o5_noopt.png", plot, AR=2/1, width=900)
displaygui(plot)

#Optimized zeros
println("\t\t\tOptimized zeros\n")
NTF5_opt = synthesizeNTF(5, OSR, opt=1)
plot = plotNTF(NTF5_opt, OSR, color=:red)
plot.title = "5th-Order Modulator (Optimized Zeros)"
saveimage(:png, "dsdemo1_o5_zopt.png", plot, AR=2/1, width=900)
displaygui(plot)

#Compare (overlay results)
plot = plotNTF(NTF5_noopt, OSR, color=:blue)
plot = plotNTF!(plot, NTF5_opt, OSR, color=:red)
plot.title = "5th-Order Modulator (Optimized Zeros - Overlay)"
saveimage(:png, "dsdemo1_o5_cmp.png", plot, AR=2/1, width=900)
displaygui(plot)


#==8th-order bandpass modulator
===============================================================================#
println("\t\tNTF Synthesis-- 8th-Order Bandpass Modulator\n")
OSR = 64
f0 = 0.125 #fs/8

NTF8bp_opt2 = synthesizeNTF(8, OSR, opt=2, f0=f0)
plot = plotNTF(NTF8bp_opt2, OSR, color=:blue)
plot.title = "8th-Order Bandpass Modulator"
saveimage(:png, "dsdemo1_o5_bp.png", plot, AR=2/1, width=900)
displaygui(plot)

:END_OF_DEMO
