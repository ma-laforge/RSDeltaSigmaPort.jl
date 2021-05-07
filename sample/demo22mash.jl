# Demonstrate MASH system
using RSDeltaSigmaPort
using RSDeltaSigmaPort.EasyPlot #set, cons
import RSDeltaSigmaPort: BoundingBox
import Printf: @sprintf
j=im


#==
===============================================================================#
ABCD = [
	1 0 0 0 1 -1  0;
	1 1 0 0 0 -2  0;
	0 1 1 0 0  0 -1;
	0 0 1 1 0  0 -2;
	0 1 0 0 0  0  0;
	0 0 0 1 0  0  0;
]
nlev = [9 9]
(ntf, stf) = calculateTF(ABCD, [1,1])
println()
@error("Need to get minreal() calcualtion to work in calculateTF() to proceed.")
println()
@info("TODO: Port rest of code once this is done.")
throw(:ERROR)
ncf1 = -ntf[2,1] #2*(z-0.5)^2/z^4
ncf2 = ntf[1,1] #(z-1)^2/z^2
#stf_eff = stf[1,1]*ncf1 + stf[2,1]*ncf2
stf_eff = cancelPZ( _zpk(_tf(stf[1,1]*ncf1) + tf(stf[2,1]*ncf2)) )
ntf_eff = ntf(2,2)*ncf2 #(z-1)^4/z^4

:END_OF_DEMO
