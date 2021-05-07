#Example lowpass/bandpass real/quadrature modulator design
include("dsexample1_fn.jl")

opt=2
#dsm = RealDSM(order=5, OSR=32, form=:CRFB, Hinf=1.5, opt=2)
dsm = RealDSM(order=5, OSR=32, form=:CRFB, Hinf=1.5, opt=opt)
results = dsexample1(dsm, LiveDemo=true)

:END_OF_EXAMPLE
