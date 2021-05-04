#RSDeltaSigmaPort: Time domain operations/utilities
#-------------------------------------------------------------------------------

"""`y = sinc_decimate(x, m, r)`

Decimate x by m-th order sinc of length r.
"""
function sinc_decimate(x, m::Int, r::Int)
	x = x[:] #Vector form
	for i in 1:m
		x = cumsum(x)
		x = vcat(x[1:r], x[r+1:end] .- x[1:end-r])/r
	end
	y = x[r:r:end]
	return y
end

#Last line
