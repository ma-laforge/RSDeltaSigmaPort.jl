#RSDeltaSigmaPort: Windowing functions
#-------------------------------------------------------------------------------

"""`w = ds_hann(n)`

A Hann window of length n. Does not smear tones located exactly in a bin.
"""
function ds_hann(n::Int)
	w = .5*(1 .- cos.(2*pi*collect(0:n-1)/n) )
	return w
end


#Last line
