# Demonstrate Demonstrate simulateDSM and simulateSNR
using RSDeltaSigmaPort
j=im

println("\t\t\tDiscrete-Time Simulation")

OSR = 32
NTF = synthesizeNTF(5, OSR, opt=1)
N = 8192
fB = ceil(Int, N/(2*OSR)); ftest=floor(Int, 2/3*fB)
u = 0.5*sin.(2π*ftest/N*collect(0:N-1)) # half-scale sine-wave input

@info("Starting simulation..."); flush(stdout); flush(stderr)
v,xn,xmax,y = simulateDSM(u, NTF)
println("\tdone.")

@info("Plotting modulator signals."); flush(stdout); flush(stderr)
plot = plotModTransient(u, v, y)
displaygui(plot)

@info("Plotting spectrum of modulator output, v."); flush(stdout); flush(stderr)
plot = plotModSpectrum(v, NTF, fB, ftest-2,
	title="Modulator Output Spectrum @ OSR = $OSR."
)
displaygui(plot)

@info("Plotting SNR vs input power."); flush(stdout); flush(stderr)
sig = v
plot = plotSNR(sig, NTF, OSR,
	title="SNR curve- theory and simulation"
)
displaygui(plot)

:END_OF_DEMO
