#Design example for a continuous-time lowpass ΔΣ ADC
using RSDeltaSigmaPort
using RSDeltaSigmaPort.EasyPlot #set, cons
import Printf: @sprintf
j=im


#==Baseband modulator (continuous-time implementation)
===============================================================================#
nlev=2
dsm = RealDSM(order=3, OSR=100, M=nlev-1, f0=0, opt=0, Hinf=1.3, form=:FB)
#dsm = RealDSM(order=3, OSR=32, M=nlev-1, f0=0, opt=0, Hinf=1.5, form=:FB)
umax = 0.83

#Parameters for the continuous-time implementation:
tdac = [0 1] #DAC timing. [0 1] means zero-delay non-return-to-zero

println("\n*** $(dsm.order)th-Order Continuous-Time Lowpass Example")

@info("Performing NTF synthesis..."); flush(stdout); flush(stderr)
#-------------------------------------------------------------------------------
NTF0 = synthesizeNTF(dsm) #Optimized zero placement
#println("p:"); display(NTF0.p)
#println("z:"); display(NTF0.z)
plotcoll_NTF = documentNTF(dsm, NTF0)
println("\tdone.")
saveimage(:png, "dsexample2_NTF.png", plotcoll_NTF, AR=2/1, width=900)
displaygui(plotcoll_NTF)
@show length(plotcoll_NTF.plotlist)
#displaygui(plotNTF(NTF0, dsm.OSR)) #ADDED MALaforge


@info("Performing time-domain simulations..."); flush(stdout); flush(stderr)
#-------------------------------------------------------------------------------
#Example spectrum
#ftest = 0.3 * 0.5/dsm.OSR/3 #~1/3 of the way across the passband
ftest = (0.5/dsm.OSR)/3 #~1/3 of the way across the passband
plot = plotExampleSpectrum(dsm, NTF0, ftest=ftest,N=2^16)
#title('Example Spectrum');
saveimage(:png, "dsexample2_PSD.png", plot, AR=2/1, width=900)
#displaygui(plot)

#SQNR plot
snrinfo = calcSNRInfo(dsm, NTF=NTF0)
plot = plotSNR(snrinfo, dsm)
	set(plot, xyaxes=set(xmin=-100, ymin=0, ymax=100))
println("\tdone.")
saveimage(:png, "dsexample2_SQNR.png", plot, AR=2/1, width=900)
#displaygui(plot)


@info("Mapping to continuous-time..."); flush(stdout); flush(stderr)
#-------------------------------------------------------------------------------
(ABCDc,tdac2) = realizeNTF_ct(NTF0, dsm.form, tdac)
(Ac, Bc, Cc, Dc) = partitionABCD(ABCDc)
display(ABCDc); map(display, [Ac, Bc, Cc, Dc])
#=MA: Seems to be a problem with ABCDc matrix:
G=0 in realizeNTF_ct() creates NaNs in ABCDc. Not sure this is expected.
=#
sys_c = _ss(Ac,Bc,Cc,Dc)
println("\tdone.")
#Verify that the sampled pulse response of the CT loop filter
#matches the impulse response of the DT prototype
n_imp = 10
y = -impL1(NTF0, n_imp) #Negate impL1 to make the plot prettier
	y_max = maximum(real.(y))
plot = plotLollipop(0:n_imp, y)
	plot.title = "Loop filter pulse/impulse responses (negated)"
yl = floor(0.9*y_max); color=:blue
	wfrm = waveform([0.5, 2], [yl, yl])
simglyph = cons(:a, glyph=set(shape=:o, size=1.5, color=color, fillcolor=color))
push!(plot,
	cons(:wfrm, wfrm, simglyph, line=set(style=:solid, color=color, width=2), label=""),
	cons(:atext, "   discrete-time", x=2, y=yl, align=:bl),
)
dt = 1/16
#yy = -pulse(sys_c, [0 0;tdac], dt, n_imp*1.0)
@warn("MA-FIXME: sys_c seems incorrectly calculated.")
t = 0:dt:n_imp
yl = floor(0.7*y_max); color=:green
	wfrm = waveform([0.5, 2], [yl, yl])
push!(plot,
	#Add yy when it works
	cons(:wfrm, wfrm, line=set(style=:solid, color=color, width=2), label=""),
	cons(:atext, "   continuous-time", x=2, y=yl, align=:bl),
)
displaygui(plot)
#throw(:STOP)

#Map the cts system to its discrete time equivalent and check the NTF
sys_d = mapCtoD(sys_c, t=tdac)
ABCD = [sys_d.a sys_d.b; sys_d.c sys_d.d]
(NTF, G) = calculateTF(ABCD)
NTF = cancelPZ(NTF)
[plot_NTF, plot_fresp] = plotcoll_NTF.plotlist
plotPZ!(plot_NTF, NTF, color=:cyan)
#Also plot the STF
LF = zpk(sys_c)
L0 = LF[1]
f = range(0, stop=0.5, N=10)
G = evalTFP(L0, NTF, f)
Gw = waveform(f, dbv.(G))
push!(plot_fresp,
	cons(:wfrm, Gw, line=set(style=:solid, color=:magenta, width=2), label=""),
)
plotcoll_NTF.title = "NTF and STF"

@info("Performing dynamic range scaling..."); flush(stdout); flush(stderr)
#-------------------------------------------------------------------------------
#!!! This code assumes that the scale factors for the DT equivalent apply 
#!!! to the CT system. A system with an RZ DAC will have inter-sample peaks 
#!!! that exceed the values an the sampling instants.
(ABCDs, umax, S) = scaleABCD(ABCD, nlev=nlev, f0=dsm.f0, xlim=1, umax=umax, N0=1e4)
S = S[1:order,1:order] #Don't worry about the extra states used in the d-t model
Sinv = inv(S)
Acs=S*Ac*Sinv; Bcs=S*Bc;  Ccs=Cc*Sinv
ABCDcs = [Acs Bcs; Ccs Dc]
sys_cs = _ss(Acs, Bcs, Ccs, Dc)
println("\tdone.")
#ABCDcs needs to be checked with CT simulations to
# 1. Verify pulse response
# 2. Verify signal swings

adc = Dict{Symbol, Any}(
	:order => dsm.order,
	:OSR => dsm.OSR,
	:opt => dsm.opt,
	:M => dsm.M,
	:f0 => dsm.f0,
	:NTF => NTF,
	:ABCD => ABCD,
	:umax => umax,
	:peak_snr => snrinfo.peak[2],
	:form => dsm.form,
	:ABCDc => ABCDc,
	:L0  => L0,
	:sys_c => sys_c,
	:ABCDcs => ABCDcs,
	:sys_cs => sys_cs,
)

:END_OF_EXAMPLE
