# Demonstrate synthesizeNTF
using RSDeltaSigmaPort
using CMDimData
using CMDimData.MDDatasets
using CMDimData.EasyPlot
CMDimData.@includepkg EasyPlotInspect
j=im

if !@isdefined(liveDemo)
	global liveDemo=0
end

fig1pos1 = [9 630 200 200]
fig1pos2 = [10 407 480 420]
fig2pos1 = [239 630 450 200]
fig2pos2 = [241 341 523 485]

println("\t\tNTF Synthesis-- 5th-order modulator\n")
OSR = 32
H = synthesizeNTF(order=5, osr=OSR, opt=0)


#Linear plot
𝑓 = vcat(range(0, stop=0.75/OSR, length=100), range(0.75/OSR, stop=0.5, length=100))
z = exp.(2j*π*𝑓)
H1 = synthesizeNTF(order=5, osr=OSR, opt=1)
magH1 = dbv(evalTF(H1, z))
magH1 = DataF1(𝑓, magH1)

#Log plot
𝑓start = 0.01
𝑓 = collect(range(𝑓start, stop=1.2, length=200)/(2*OSR))
z = exp.(2j*π*𝑓)
𝑓norm = 𝑓*(2*OSR)
magH1 = dbv(evalTF(H1, z))
magH1 = DataF1(𝑓norm, magH1)

#==Generate EasyPlot
===============================================================================#
alabels = cons(:a, labels=set(xaxis="Normalized frequency (1→f_s)", yaxis="dB"))
linlin = cons(:a, xyaxes=set(xscale=:lin, yscale=:lin, ymin=-100)) #Limit ymin to finite value
loglin = cons(:a, xyaxes=set(xscale=:log, yscale=:lin, xmin=10^-2, ymin=-100)) #Limit ymin to finite value

plot = push!(cons(:plot, loglin, alabels, title = "NTF Magnitude Response", legend=false),
	cons(:wfrm, magH1, label="|H1|"),
)

pcoll = push!(cons(:plot_collection, title="Sample Plot"), plot)
EasyPlot.displaygui(:InspectDR, pcoll)

:END_OF_DEMO
