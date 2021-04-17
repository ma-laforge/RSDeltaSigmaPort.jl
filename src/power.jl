#RSDeltaSigmaPort: Compute power values.
#-------------------------------------------------------------------------------

"""`dbm(v,R=50) = 10*log10(v^2/R*1000)`

The equivalent in dBm of an RMS voltage v
"""
dbm(v; R::Float64=50.0) = 10*log10.(abs.(v .^ 2)/R) .+ 30;

"""`dbp(x) = 10*log10(x)`

The dB equivalent of the power x
"""
dbp(x) = 10*log10.(abs.(x))

"""`dbv(x) = 20*log10(abs(x))`

The dB equivalent of the voltage x
"""
dbv(x) = 20*log10.(abs.(x))

"`y = rms(x; no_dc=false)`"
function rms(x; no_dc::Bool=false)
	if no_dc
		x = x - mean(x)
	end
	return norm(x)/sqrt(length(x))
end

"""`v = undbm(p, z=50) = sqrt(z*10^(p/10-3))`

RMS voltage equivalent to a power p in dBm
"""
undbm(p; z::Float64=50.0) = sqrt.(z*10 .^ (p/10 .- 3))

"""`y = undbp(x)`

Convert x from dB to a power
"""
undbp(x) = 10 .^ (x/10)

"""`y = undbv(x)`

Convert x from dB to a voltage
"""
undbv(x) = 10 .^ (x/20)

#Last line
