# Demonstrate synthesizeNTF
using RSDeltaSigmaPort
using CMDimData
using CMDimData.EasyPlot
CMDimData.@includepkg EasyPlotInspect
j=im

#==5th-order modulator
===============================================================================#
println("\t\tNTF Synthesis-- 5th-Order Modulator\n")

OSR = 32

NTF = synthesizeNTF(5, OSR, opt=0)
pcoll = plotNTF(NTF, OSR, color=:blue)
pcoll.title = "5th-Order Modulator"
plot_gui = EasyPlot.displaygui(:InspectDR, pcoll)

println("\t\t\tOptimized zeros\n")
NTF = synthesizeNTF(5, OSR, opt=1)
pcoll = plotNTF(NTF, OSR, color=:red)
pcoll.title = "5th-Order Modulator (Optimized Zeros)"
plot_gui = EasyPlot.displaygui(:InspectDR, pcoll)

NTF = synthesizeNTF(5, OSR, opt=0)
pcoll = plotNTF(NTF, OSR, color=:blue)
NTF = synthesizeNTF(5, OSR, opt=1)
pcoll = plotNTF!(pcoll, NTF, OSR, color=:red)
pcoll.title = "5th-Order Modulator (Optimized Zeros - Overlay)"

plot_gui = EasyPlot.displaygui(:InspectDR, pcoll)
ploth = 800; plotw = round(Int, ploth*2)
EasyPlot._write(:png, "dsdemo1_1.png", plot_gui, dim=set(w=plotw, h=ploth))


#==8th-order bandpass modulator
===============================================================================#
OSR = 64
f0 = 0.125 #fs/8
println("\t\tNTF Synthesis-- 8th-Order Bandpass Modulator\n")
NTF = synthesizeNTF(8, OSR, opt=2, f0=f0)
pcoll = plotNTF(NTF, OSR, color=:blue)
pcoll.title = "8th-Order Bandpass Modulator"
plot_gui = EasyPlot.displaygui(:InspectDR, pcoll)
ploth = 800; plotw = round(Int, ploth*2)
EasyPlot._write(:png, "dsdemo1_bp.png", plot_gui, dim=set(w=plotw, h=ploth))

:END_OF_DEMO
