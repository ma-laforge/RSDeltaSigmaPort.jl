#RSDeltaSigmaPort: Windowing functions
#-------------------------------------------------------------------------------

"""`w = ds_hann(n)`

A Hann window of length n. Does not smear tones located exactly in a bin.
"""
function ds_hann(n::Int)
	w = .5*(1 .- cos.(2*pi*collect(0:n-1)/n) )
	return w
end


"""`circ_smooth(x,n)`
"""
function circ_smooth(x::Array{T, 2}, n::Int) where T <: Number
	szx = size(x)
	if szx[1] != 1
		throw("Arrays with nrows!=1 are not yet supported")
	end
	nx = szx[2]
	w = ds_hann(n)'/(n/2)
	xw = conv(x, w)
	y = circshift(hcat(xw[:, n:nx], xw[:, 1:n-1] .+ xw[:, nx+1:end]), (0, div(n,2)-1))
	return y
end

function circ_smooth(x::Vector, n::Int)
	x = conv2seriesmatrix2D(x)
	return circ_smooth(x, n)[:] #Keep things in vector form
end


#Last line
