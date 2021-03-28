# Demonstrate Demonstrate simulateDSM and simulateSNR
using RSDeltaSigmaPort
j=im

println("\t\t\tDiscrete-Time Simulation")

OSR = 32
H = synthesizeNTF(5, OSR, opt=1)
N = 8192
fB = ceil(N/(2*OSR)); ftest=floor(2/3*fB)
u = 0.5*sin.(2π*ftest/N*collect(0:N-1)) # half-scale sine-wave input
throw(:NOT_YET_IMPLEMENTED)
v = simulateDSM(u,H)


:END_OF_DEMO
