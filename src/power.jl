#RSDeltaSigmaPort: Compute power values.
#-------------------------------------------------------------------------------

"""`dbm(v,R=50) = 10*log10(v^2/R*1000)`
The equivalent in dBm of an rms voltage v
"""
dbm(v; R::Float64=50.0) = 10*log10.(abs.(v .^ 2)/R) .+ 30;

"""`dbp(x) = 10*log10(x)`
the dB equivalent of the power x
"""
dbp(x) = 10*log10.(abs.(x))

"""`dbv(x) = 20*log10(abs(x))`
the dB equivalent of the voltage x
"""
dbv(x) = 20*log10.(abs.(x))

"`y = rms(x; no_dc=false)`"
function rms(x; no_dc::Bool=false)
	if no_dc
		x = x - mean(x)
	end
	return norm(x)/sqrt(length(x))
end

#Last line
