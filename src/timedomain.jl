#RSDeltaSigmaPort: Time domain operations/utilities
#-------------------------------------------------------------------------------

"""`y = delay(x, n=1)`

Delay signal x by n samples
"""
function delay(x, n::Int=1)
	local y
	nx = length(x)
	if nx <= n
		y = zeros(size(x))
	else
		if 1 == length(size(x))
			y = vcat(zeros(n), x[1:nx-n])
		elseif 1 == size(x, 2) #Still a row vector, but 2D+ now
			y = vcat(zeros(n,1), x[1:nx-n,[1]])
		else #Assume a single row vector:
			y = hcat(zeros(1,n), x[[1],1:nx-n])
		end
	end
	return y
end


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

"""`applywnd(v, wnd)`

Safely apply window `wnd`.

Avoids calling `fft()` on an accidentally large matrix.
"""
function applywnd(v, wnd)
	szv = size(v); szw = size(wnd)

	if length(szv) < 1 || length(szw) < 1
		throw("applywnd: Does not support 0-length arrays")
	elseif length(szv) > 2 || length(szw) > 2
		throw("applywnd: Does not yet support 3D+ arrays")
	elseif szv == szw
		return v .* wnd
	elseif 1 == length(szw) #Wnd is 1D (column) matrix
		if szv[2] == szw[1]
			return v .* wnd'
		end
	elseif 1 == szw[1] #Wnd is a row matrix
		if szv[2] == szw[2]
			return v .* wnd
		end
	end

	throw(ArgumentError("Unexpected input sizes: $szv & $szw"))
end

#Last line
