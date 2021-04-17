#RSDeltaSigmaPort: Generate time-domain plots of ΔΣ signals
#-------------------------------------------------------------------------------


#==Waveform builders
===============================================================================#
function _stairs__(x::Vector, y::Vector, sweepid::String)
	validatexy_vec(x, y)
	if length(x) < 1
		throw("Does not support zero-length arrays")
	end
	xs = similar(x, 2*length(x))
	ys = similar(y, 2*length(y))

	ii = io = 1
	while ii < length(x)
		xs[io] = x[ii]
		xs[io+1] = x[ii+1]
		ii += 1; io +=2
	end
	if ii <= length(x) #last point
		Δx = 1 #Step size (if x has a single point)
		if length(x) > 1; Δx = x[end]-x[end-1]; end
		xs[io] = x[ii]
		xs[io+1] = x[ii]+Δx #Use final step size
	end

	ii = io = 1
	while ii <= length(y)
		_y = y[ii]
		ys[io] = _y
		ys[io+1] = _y
		ii += 1; io +=2
	end
	return DataF1(xs, ys)
end

function _stairs__(x::Array{Float64,2}, y::Array{Float64,2}, sweepid::String)
	validatexy_1param(x, y)
	N = size(x, 1)
	wfrm = fill(DataRS{DataF1}, PSweep(sweepid, collect(1:N))) do i
		return _stairs__(x[i,:], y[i,:], "") #Returns from implicit "do" function
	end
	return wfrm
end

"""`wfrm_stairs(x, y; sweepid="i")`

Create a staircase-waveform object with a given string to identify the sweep.
"""
wfrm_stairs(x, y; sweepid::String="i") = _stairs__(x, y, sweepid)


#==Time domain plots of modulator
===============================================================================#
function plotModTransient(inputSig, outputSig, otherSig...; legend::Bool=true, color=:blue)
	plot = cons(:plot, linlin, title = "Modulator Input & Output", legend=legend,
		labels=set(xaxis="Sample Number", yaxis="Amplitude [V]"),
	)

	#Ensure we have 2D-Array{}s, not Vector{}s:
	inputSig = conv2seriesmatrix2D(inputSig)
	outputSig = conv2seriesmatrix2D(outputSig)

	#Generate matrix of sample numbers:
	M = size(inputSig,1) #Number of rows (ex: quantizers)
	N = size(inputSig,2) #Number of samples
	n = 1:N
	sn = ones(Float64, M)*n' #Matrix of "sample number" (Use Float64s for downstream functions)

	inputSig = wfrm_stairs(sn, inputSig)
	outputSig = wfrm_stairs(sn, outputSig)

	push!(plot,
		cons(:wfrm, outputSig, line=set(style=:solid, color=color, width=2), label="output"),
		cons(:wfrm, inputSig, line=set(style=:solid, color=:red, width=2), label="input"),
#		cons(:wfrm, otherSig, line=set(style=:solid, color=:black, width=2), label="y"),
	)
	return plot
end

#Last line
