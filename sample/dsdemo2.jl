# Demonstrate simulateDSM, (simulateSNR and predictSNR) => plotSNR
using RSDeltaSigmaPort
using RSDeltaSigmaPort.EasyPlot #set, cons
import RSDeltaSigmaPort: BoundingBox
import Printf: @sprintf
j=im


#==Baseband modulator
===============================================================================#
println("\n*** 5th order, 2-level, baseband modulator")
OSR = 32
N = 8192
dsm = RealDSM(order=5, OSR=OSR, opt=1)

@info("Performing ΔΣ simulation..."); flush(stdout); flush(stderr)
#-------------------------------------------------------------------------------
ftest = dsm.f0 + 1/(3*OSR); fband = default_fband(dsm)
fband[1] += 2/N #Ignore first 2 low-frequency points
(u, iftest) = genTestTone(dsm, dbv(0.5), ftest, N=N)
NTF = synthesizeNTF(dsm)
simresult = simulateDSM(u, NTF)
println("\tdone.")

@info("Plotting modulator signals")
#-------------------------------------------------------------------------------
plot = plotModTransient(u, simresult.v)
	set(plot, xyaxes=set(xmin=0, xmax=300, ymin=-1.2, ymax=1.2))
saveimage(:png, "dsdemo2_o5_sig.png", plot, AR=2/1, width=900)
displaygui(plot)

@info("Plotting output spectrum (simulated vs theory)")
#-------------------------------------------------------------------------------
specinfo = calcSpecInfo(simresult.v, NTF, fband, ftest)
plot = plotModSpectrum(specinfo)
	plot.title="Modulator Output Spectrum @ OSR = $OSR."
saveimage(:png, "dsdemo2_o5_spec.png", plot, AR=2/1, width=900)
displaygui(plot)

@info("Plotting SNR vs input power")
#-------------------------------------------------------------------------------
plot = plotSNR(simresult.v, NTF, OSR,
	title="SNR curve- theory and simulation"
)
saveimage(:png, "dsdemo2_o5_snr.png", plot, AR=2/1, width=900)
displaygui(plot)


#==Bandpass modulator
===============================================================================#
println("\n*** 8th order, 2-level, bandpass modulator")
OSR=64
f0 = 1/8
N = 8192
dsm = RealDSM(order=8, OSR=OSR, f0=f0, opt=1)
display(dsm)

@info("Performing ΔΣ simulation..."); flush(stdout); flush(stderr)
#-------------------------------------------------------------------------------
ftest = dsm.f0 + 1/(6*OSR); fband = default_fband(dsm)
(u, iftest) = genTestTone(dsm, dbv(0.5), ftest, N=N) #Half-scale sine-wave input
NTF = synthesizeNTF(dsm)
simresult = simulateDSM(u, NTF)
println("\tdone.")

@info("Plotting modulator signals")
#-------------------------------------------------------------------------------
plot = plotModTransient(u, simresult.v)
	set(plot, xyaxes=set(xmin=0, xmax=300, ymin=-1.2, ymax=1.2))
saveimage(:png, "dsdemo2_bp_o5_sig.png", plot, AR=2/1, width=900)
displaygui(plot)

@info("Plotting output spectrum (simulated vs theory)")
#-------------------------------------------------------------------------------
specinfo = calcSpecInfo(simresult.v, NTF, fband, ftest)
plot = plotModSpectrum(specinfo)
	plot.title="Modulator Output Spectrum @ OSR = $OSR."
saveimage(:png, "dsdemo2_bp_o5_spec.png", plot, AR=2/1, width=900)
displaygui(plot)

@info("Plotting SNR vs input power")
#-------------------------------------------------------------------------------
plot = plotSNR(simresult.v, NTF, OSR, f0=f0,
	title="SNR curve- theory and simulation"
)
saveimage(:png, "dsdemo2_bp_o5_snr.png", plot, AR=2/1, width=900)
displaygui(plot)


#==15-step baseband modulator
===============================================================================#
println("\n*** 7th order, 15-step, baseband modulator")
OSR = 8
M = 16 #Shouldn'nt this be 15 (nlev=16?)
N = 8192
Hinf_list = [2.0, 8.0]
dsm_list = [RealDSM(order=7, OSR=OSR, M=M, opt=1, Hinf=Hinf) for Hinf in Hinf_list]
color_list = [:blue, :green]

@info("Performing ΔΣ simulation for H(∞)=$(Hinf_list)..."); flush(stdout); flush(stderr)
#-------------------------------------------------------------------------------
id_list = [@sprintf("H(∞)=%.1f", Hinf) for Hinf in Hinf_list]
ftest = 1/(7*OSR); fband = default_fband(OSR)
fband[1] += 2/N #Ignore first 2 low-frequency points
(u, iftest) = genTestTone(dsm, dbv(0.5*M), ftest, N=N) #Half-scale sine-wave input
NTF_list = [synthesizeNTF(dsm) for dsm in dsm_list]
v_list = [simulateDSM(u, NTF, nlev=M+1).v for NTF in NTF_list] #simulateDSM(..).v is mod output, v.
println("\tdone.")

@info("Plotting input/output characteristics of ΔΣ simulations")
#-------------------------------------------------------------------------------
ioplotc = cons(:plot_collection, title="15-step / 7th-order / H(∞)=$(Hinf_list)")

#Plot input & output transients:
for (v, c) in zip(v_list, color_list) #Each simulated output
	local plot = plotModTransient(u, v, color=c, legend=false)
		set(plot, xyaxes=set(xmin=0, xmax=100, ymin=-16, ymax=16))
	push!(ioplotc, plot)
end

#Append SNR vs input curves:
plot = plotSNR(v_list[1], NTF_list[1], OSR, nlev=M+1, color=color_list[1], legend=false, title="SQNR")
plotSNR!(plot, v_list[2], NTF_list[2], OSR, nlev=M+1, color=color_list[2])
	set(plot, xyaxes=set(xmin=-100, xmax=0, ymin=0, ymax=120))
push!(ioplotc, plot)

#Specify plot locations to help readability:
ioplotc.bblist = [ #Format plot locations
	BoundingBox(0, 0.5, 0, 0.5), #1st modulator transient
	BoundingBox(0, 0.5, 0.5, 1), #2nd modulator transient
	BoundingBox(0.5, 1, 0, 1), #SNR curve
]
saveimage(:png, "dsdemo2_o7_15s_io.png", ioplotc, AR=2/1, width=900)
displaygui(ioplotc)

@info("Plotting output spectrum (simulated vs theory)")
#-------------------------------------------------------------------------------
plot = plotModSpectrum()
	plot.title="Modulator Output Spectrum @ OSR = $OSR."
	set(plot, xyaxes=set(ymin=-160))
for i in keys(NTF_list)
	local specinfo = calcSpecInfo(v_list[i], NTF_list[i], fband, ftest, M=M)
	plotModSpectrum!(plot, specinfo, id=id_list[i], color=color_list[i])
end
saveimage(:png, "dsdemo2_o7_15s_spec.png", plot, AR=2/1, width=900)
displaygui(plot)

:END_OF_DEMO
