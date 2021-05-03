#RSDeltaSigmaPort: Base plot generation tools
#-------------------------------------------------------------------------------


#==Useful constants
===============================================================================#
const linlin = cons(:a, xyaxes=set(xscale=:lin, yscale=:lin))
const linlog = cons(:a, xyaxes=set(xscale=:lin, yscale=:log))
const loglin = cons(:a, xyaxes=set(xscale=:log, yscale=:lin))
const loglog = cons(:a, xyaxes=set(xscale=:log, yscale=:log))

const dfltline = cons(:a, line=set(style=:solid, color=:blue, width=2))


#==Waveform builders
===============================================================================#
function _wfrm__(x::Vector, y::Vector, sweepid::String)
	validatexy_vec(x, y)
	return DataF1(x,y)
end
_wfrm__(x::AbstractVector, y::Vector, sweepid::String) = _wfrm__(collect(x), y, sweepid)
function _wfrm__(x::Array{Float64,2}, y::Array{Float64,2}, sweepid::String)
	validatexy_1param(x, y)
	N = size(x, 1)
	wfrm = fill(DataRS{DataF1}, PSweep(sweepid, collect(1:N))) do i
		return DataF1(x[i,:], y[i,:]) #Returns from implicit "do" function
	end
	return wfrm
end

"""`waveform(x, y; sweepid="i")`

Create a waveform object with a given string to identify the sweep.
"""
waveform(x, y; sweepid::String="i") = _wfrm__(x, y, sweepid)


#==logsmooth
===============================================================================#
"""`(f, p) = logsmooth(X,inBin,nbin=8,n=3)`

Smooth the fft, X, and convert it to dB.

 - Use nbin bins from 0 to 3*inBin, thereafter increase bin sizes by a factor
   of 1.1, staying less than 2^10.
 - For the n sets of bins inBin+[0:2], 2*inBin+[0:2], ... n*inBin+[0:2], don't
   do averaging. This way, the noise BW and the scaling of the tone and its
   harmonics are unchanged.

Unfortunately, harmonics above the nth appear smaller than they 
really are because their energy is averaged over many bins.
"""
function logsmooth(X, inBin; nbin::Int=8, n::Int=3)
	N=length(X); N2=div(N,2)
	f1 = mod(inBin-1, nbin) + 1
	startbin = vcat(f1:nbin:inBin-1, inBin:inBin+2)
	for i in 1:n
		startbin = vcat(startbin, startbin[end]+1:nbin:i*inBin-1, i*inBin .+ (0:2))
	end
	m = startbin[length(startbin)]+nbin
	while m < N2
		startbin = vcat(startbin, m)
		nbin = min(nbin*1.1, 2^10)
		m = round(Int, m+nbin)
	end
	stopbin = vcat(startbin[2:length(startbin)] .- 1, N2)
	f = ((startbin .+ stopbin)/2 .- 1)/N
	p = zeros(size(f))

	for i in 1:length(f)
		p[i] = dbp(norm(X[startbin[i]:stopbin[i]])^2 / (stopbin[i]-startbin[i]+1))
	end

	return (f, p)
end

#Last line
