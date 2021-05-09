#Demonstrate the designLCBP function
#=WARNING:
designLCBP is very ill-behaved due to deficiencies in MATLAB's constr() and
minimax(). These functions frequently take a viable set of parameters as a
starting point and turn them into an unstable modulator. I should provide a
more robust implementation...
=#
using RSDeltaSigmaPort
using RSDeltaSigmaPort.EasyPlot #set, cons
import Printf: @sprintf
j=im


#==
===============================================================================#
println("\n*** Continuous-Time LC Modulator Design")

n = 3
OSR = 64
opt = 2
Hinf = 1.7
f0 = 1/16
t = [0.5 1]
form = :FB
dbg = true
(param, H, L0, ABCD) = designLCBP(n, OSR=OSR, opt=opt,
	Hinf=Hinf, f0=f0, t=t, form=form, dbg=dbg
)

:END_OF_DEMO
