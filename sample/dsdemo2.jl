# Demonstrate Demonstrate simulateDSM and simulateSNR
using RSDeltaSigmaPort
j=im


#==Discrete-Time Simulation
===============================================================================#
println("\n\t\t\tDiscrete-Time Simulation")

OSR = 32
N = 8192
fB = ceil(Int, N/(2*OSR)); ftest=floor(Int, 2/3*fB)
u = 0.5*sin.(2π*ftest/N * (0:N-1)) # half-scale sine-wave input
NTF = synthesizeNTF(5, OSR, opt=1)

@info("Performing ΔΣ simulation..."); flush(stdout); flush(stderr)
v,xn,xmax,y = simulateDSM(u, NTF)
println("\tdone.")

@info("Plotting modulator signals."); flush(stdout); flush(stderr)
plot = plotModTransient(u, v, y)
displaygui(plot)

@info("Plotting spectrum of modulator output, v."); flush(stdout); flush(stderr)
plot = plotModSpectrum(v, NTF, 3:fB+1, ftest-2,
	title="Modulator Output Spectrum @ OSR = $OSR."
)
displaygui(plot)

@info("Plotting SNR vs input power."); flush(stdout); flush(stderr)
plot = plotSNR(v, NTF, OSR,
	title="SNR curve- theory and simulation"
)
displaygui(plot)


#==Bandpass Modulator
===============================================================================#
println("\n\t\t\tBandpass Modulator")

OSR=64
f0 = 1/8
N = 8192
fB = ceil(Int, N/(2*OSR)); ftest=round(Int, f0*N + 1/3*fB)
u = 0.5*sin.(2π*ftest/N * (0:N-1)) #half-scale sine-wave input
NTF = synthesizeNTF(8, OSR, opt=1, f0=f0)

@info("Performing ΔΣ simulation..."); flush(stdout); flush(stderr)
v,xn,xmax,y = simulateDSM(u, NTF)
println("\tdone.")

@info("Plotting modulator signals."); flush(stdout); flush(stderr)
plot = plotModTransient(u, v, y)
displaygui(plot)

@info("Plotting spectrum of modulator output, v."); flush(stdout); flush(stderr)
f1 = round(Int, (f0-0.25/OSR)*N)
f2 = round(Int, (f0+0.25/OSR)*N)
plot = plotModSpectrum(v, NTF, f1:f2, ftest-f1+1,
	title="Modulator Output Spectrum @ OSR = $OSR."
)
displaygui(plot)

@info("Plotting SNR vs input power."); flush(stdout); flush(stderr)
plot = plotSNR(v, NTF, OSR, f0=f0,
	title="SNR curve- theory and simulation"
)
displaygui(plot)

:END_OF_DEMO
