# Realization and dynamic range scaling
using RSDeltaSigmaPort
using RSDeltaSigmaPort.EasyPlot #set, cons
import Printf: @sprintf
j=im

function stringify(vA, digits=6)
	#ctx=IOContext(stdout, :compact=>true)
	#IOContext does not support # of digits
	#_s(v) = @sprintf("%.6f", v) #Cannot parameterize format specifier
	return [round(v, digits=6) for v in vA]
end


#==
===============================================================================#
println("\n*** Modulator realization and scaling")
OSR=42
order=5; opt=1

NTF = synthesizeNTF(order, OSR, opt=opt)
plot = plotNTF(NTF, OSR)
plot.title = string(order, "th-Order Modulator")
#saveimage(:png, "dsdemo3_ntf.png", plot, AR=2/1, width=900)
displaygui(plot)

@info("Realizing NTF..."); flush(stdout); flush(stderr)
#-------------------------------------------------------------------------------
a,g,b,c = realizeNTF(NTF)
b = [b[1] zeros(1,length(b)-1)] #Use a single feed-in for the input

println("\nUnscaled modulator")
println("   DAC feedback coefficients = ", join(stringify(a), ", "))
println("   resonator feedback coefficients = ", join(stringify(g), ", "))

println()
@info("Calculating the state maxima...")
#-------------------------------------------------------------------------------
ABCD = stuffABCD(a,g,b,c)
u = range(0, stop=0.6, length=30)
N = 10^4
plot = plotStateMaxima(u, ABCD, N=N)
displaygui(plot)

println()
@info("Calculating the scaled coefficients...")

#-------------------------------------------------------------------------------
ABCDs, umax = scaleABCD(ABCD, N=N)
as, gs, bs, cs = mapABCD(ABCDs)
println("\nScaled modulator, umax=", @sprintf("%.2f", umax))
println("   DAC feedback coefficients = ", join(stringify(as), ", "))
println("   resonator feedback coefficients = ", join(stringify(gs), ", "))
println("   interstage coefficients = ", join(stringify(cs), ", "))
println("   feed-in coefficients = ", join(stringify(bs), ", "))

println()
@info("Calculating the state maxima...")
#-------------------------------------------------------------------------------
u = range(0, stop=umax, length=30)
N = 10^4

plot = plotStateMaxima(u, ABCDs, N=N)
displaygui(plot)

:END_OF_DEMO
