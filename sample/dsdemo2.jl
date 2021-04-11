# Demonstrate Demonstrate simulateDSM and simulateSNR
using RSDeltaSigmaPort
j=im

println("\t\t\tDiscrete-Time Simulation")

OSR = 32
NTF = synthesizeNTF(5, OSR, opt=1)
N = 8192
fB = ceil(N/(2*OSR)); ftest=floor(2/3*fB)
u = 0.5*sin.(2π*ftest/N*collect(0:N-1)) # half-scale sine-wave input

@info("Starting simulation..."); flush(stdout); flush(stderr)
v,xn,xmax,y = simulateDSM(u, NTF)
println("\ndone.")

@info("Ploting modulator signals..."); flush(stdout); flush(stderr)
plot = plotModTransient(u, v, y)
displaygui(plot)

@info("Plotting spectrum of modulator output, v..."); flush(stdout); flush(stderr)
plot = RSDeltaSigmaPort.plotModSpectrum(v, NTF, Int(fB), Int(ftest-2),
	title="Modulator Output Spectrum @ OSR = $OSR."
)
displaygui(plot)


:END_OF_DEMO
