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
ntf0 = synthesizeNTF(dsm) #Optimized zero placement
#println("p:"); display(ntf0.p)
#println("z:"); display(ntf0.z)
plot = documentNTF(dsm, ntf0)
println("\tdone.")
saveimage(:png, "dsexample2_NTF.png", plot, AR=2/1, width=900)
displaygui(plot)

@info("Performing time-domain simulations..."); flush(stdout); flush(stderr)
#-------------------------------------------------------------------------------
#Example spectrum
#ftest = 0.3 * 0.5/dsm.OSR/3 #~1/3 of the way across the passband
ftest = (0.5/dsm.OSR)/3 #~1/3 of the way across the passband
plot = plotExampleSpectrum(dsm, ntf0, ftest=ftest,N=2^16)
#title('Example Spectrum');
saveimage(:png, "dsexample2_PSD.png", plot, AR=2/1, width=900)
displaygui(plot)

#SQNR plot
snrinfo = calcSNRInfo(dsm, NTF=ntf0)
plot = plotSNR(snrinfo, dsm)
	set(plot, xyaxes=set(xmin=-100, ymin=0, ymax=100))
println("\tdone.")
saveimage(:png, "dsexample2_SQNR.png", plot, AR=2/1, width=900)
displaygui(plot)


@info("Mapping to continuous-time..."); flush(stdout); flush(stderr)
#-------------------------------------------------------------------------------
(ABCDc,tdac2) = realizeNTF_ct(ntf0, dsm.form, tdac)
(Ac, Bc, Cc, Dc) = partitionABCD(ABCDc)
sys_c = RSDeltaSigmaPort._ss(Ac,Bc,Cc,Dc)
println("\tdone.")
#Verify that the sampled pulse response of the CT loop filter
#matches the impulse response of the DT prototype
n_imp = 10
y = -impL1(ntf0, n_imp) #Negate impL1 to make the plot prettier
plot = plotLollipop(0:n_imp, y)
displaygui(plot)


:END_OF_EXAMPLE
