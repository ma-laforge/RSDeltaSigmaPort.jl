#Lowpass and Bandpass Demonstration Script
include("dsexample1_fn.jl")
#Reminder: nlev=M+1

function wait_for_user(dsm)
	println("\nSimulation complete for modulator:\n\t$dsm")
	println("\n*** Press ENTER to continue ***")
	readline(stdin)
end

#2nd-order lowpass
dsm = RealDSM(order=2, OSR=16, M=8, Hinf=2.0, opt=0)
dsexample1(dsm, LiveDemo=true)
wait_for_user(dsm)

#5th-order lowpass
dsm = RealDSM(order=5, OSR=16, M=8, Hinf=2.0, opt=0)
dsexample1(dsm, LiveDemo=true)
wait_for_user(dsm)

#5th-order lowpass with optimized zeros
dsm = RealDSM(order=5, OSR=16, M=8, Hinf=2.0, opt=1)
dsexample1(dsm, LiveDemo=true)
wait_for_user(dsm)

#5th-order lowpass with optimized zeros and larger Hinf
dsm = RealDSM(order=5, OSR=16, M=8, Hinf=3.0, opt=1)
dsexample1(dsm, LiveDemo=true)
wait_for_user(dsm)

#7th-order lowpass; Hinf=2
dsm = RealDSM(order=7, OSR=8, M=16, Hinf=2.0, opt=1)
dsexample1(dsm, ampdB=-6.0, LiveDemo=true)
wait_for_user(dsm)

#7th-order lowpass; Hinf=8
dsm = RealDSM(order=7, OSR=8, M=16, Hinf=8.0, opt=1)
dsexample1(dsm, ampdB=-6.0, LiveDemo=true)
wait_for_user(dsm)

#6th-order bandpass
f0 = 1/6 #Normalized center frequency
dsm = RealDSM(order=6, OSR=16, M=8, f0=f0, Hinf=2.0, opt=1)
dsexample1(dsm, ampdB=-6.0, LiveDemo=true)
wait_for_user(dsm)

:END_OF_DEMO
