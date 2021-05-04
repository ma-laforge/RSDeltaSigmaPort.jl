#Design example for a continuous-time lowpass ΔΣ ADC
using RSDeltaSigmaPort
using RSDeltaSigmaPort.EasyPlot #set, cons
import Printf: @sprintf
j=im


#==Baseband modulator (continuous-time implementation)
===============================================================================#
nlev=2
#dsm = RealDSM(order=3, OSR=100, M=nlev-1, f0=0, opt=0, Hinf=1.3, form=:FB)
#dsm = RealDSM(order=3, OSR=32, M=nlev-1, f0=0, opt=0, Hinf=1.5, form=:FB)
dsm = RealDSM(order=3, OSR=32, M=nlev-1, f0=0, opt=2, Hinf=1.5, form=:FB)
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
displaygui(plotNTF(NTF0, dsm.OSR)) #ADDED MALaforge


@info("Performing time-domain simulations..."); flush(stdout); flush(stderr)
#-------------------------------------------------------------------------------
#Example spectrum
#ftest = 0.3 * 0.5/dsm.OSR/3 #~1/3 of the way across the passband
ftest = (0.5/dsm.OSR)/3 #~1/3 of the way across the passband
plot = plotExampleSpectrum(dsm, NTF0, ftest=ftest,N=2^12) #16
#title('Example Spectrum');
saveimage(:png, "dsexample2_PSD.png", plot, AR=2/1, width=900)
displaygui(plot)

#SQNR plot
snrinfo = calcSNRInfo(dsm, NTF=NTF0)
plot = plotSNR(snrinfo, dsm)
	set(plot, xyaxes=set(xmin=-100, ymin=0, ymax=100))
println("\tdone.")
saveimage(:png, "dsexample2_SQNR.png", plot, AR=2/1, width=900)
displaygui(plot)


@info("Mapping to continuous-time..."); flush(stdout); flush(stderr)
#-------------------------------------------------------------------------------
(ABCDc,tdac2) = realizeNTF_ct(NTF0, dsm.form, tdac)
(Ac, Bc, Cc, Dc) = partitionABCD(ABCDc)
sys_c = _ss(Ac,Bc,Cc,Dc)
println("\tdone.")
#Verify that the sampled pulse response of the CT loop filter
#matches the impulse response of the DT prototype
n_imp = 10
y = -impL1(NTF0, n_imp) #Negate impL1 to make the plot prettier
	y_max = maximum(real.(y))
plot = plotLollipop(0:n_imp, y, color=:blue)
	plot.title = "Loop filter pulse/impulse responses (negated)"
yl = floor(0.9*y_max); color=:blue
	ylw = waveform([0.5, 2], [yl, yl])
simglyph = cons(:a, glyph=set(shape=:o, size=1.5, color=color, fillcolor=color))
push!(plot,
	cons(:wfrm, ylw, simglyph, line=set(style=:solid, color=color, width=2), label=""),
	cons(:atext, "   discrete-time", x=2, y=yl, align=:bl),
)
dt = 1/16
t = collect(0:dt:n_imp)
yy = -pulse(sys_c, [0 0;tdac], dt, n_imp*1.0)
	yy = squeeze(yy)
	yyw = waveform(t, yy)
yl = floor(0.7*y_max); color=:green
	ylw = waveform([0.5, 2], [yl, yl])
push!(plot,
	#Add yy when it works
	cons(:wfrm, yyw, line=set(style=:solid, color=color, width=2), label=""),
	cons(:wfrm, ylw, line=set(style=:solid, color=color, width=2), label=""),
	cons(:atext, "   continuous-time", x=2, y=yl, align=:bl),
)
displaygui(plot)

#Map the cts system to its discrete time equivalent and check the NTF
(sys_d, Gp) = mapCtoD(sys_c, t=tdac)
ABCD = [sys_d.A sys_d.B; sys_d.C sys_d.D]
(NTF, G) = calculateTF(ABCD)
NTF = cancelPZ(NTF)
(plot_NTF, plot_fresp) = plotcoll_NTF.plotlist
plotPZ!(plot_NTF, NTF, color=:cyan)
#Also plot the STF
L0 = _zpk(sys_c[1,1])
@warn("FIXME: Again, numerical error in ss->zpk gives 0 gain")
f = range(0, stop=0.5, length=10)
G = evalTFP(L0, NTF, f)
Gw = waveform(f, dbv.(G))
push!(plot_fresp,
	cons(:wfrm, Gw, line=set(style=:solid, color=:magenta, width=2), label=""),
)
plotcoll_NTF.title = "NTF and STF"
displaygui(plotcoll_NTF)

@info("Performing dynamic range scaling..."); flush(stdout); flush(stderr)
#-------------------------------------------------------------------------------
#!!! This code assumes that the scale factors for the DT equivalent apply 
#!!! to the CT system. A system with an RZ DAC will have inter-sample peaks 
#!!! that exceed the values an the sampling instants.
(ABCDs, umax, S) = scaleABCD(ABCD, nlev=nlev, f=dsm.f0, xlim=1, umax=umax, N0=10^4)
S = S[1:dsm.order,1:dsm.order] #Don't worry about the extra states used in the d-t model
println("\tdone.")
println("\nScaled ABCD matrix (ABCDs):")
@show umax
display(ABCDs)

#Compute ABCDcs:
Sinv = inv(S)
Acs=S*Ac*Sinv; Bcs=S*Bc;  Ccs=Cc*Sinv
ABCDcs = [Acs Bcs; Ccs Dc]
sys_cs = _ss(Acs, Bcs, Ccs, Dc)
#ABCDcs needs to be checked with CT simulations to
# 1. Verify pulse response
# 2. Verify signal swings
println("\nABCDcs matrix:")
display(ABCDcs)

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

println()
display(adc)

:END_OF_EXAMPLE
