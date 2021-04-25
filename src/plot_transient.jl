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
function plotModTransient(; legend::Bool=true)
	plot = cons(:plot, linlin, title = "Modulator Input & Output", legend=legend,
		labels=set(xaxis="Sample Number", yaxis="Amplitude [V]"),
	)
	return plot
end

function plotModTransient!(plot, sig; color=:blue, label::String="")
	#Ensure we have a 2D-`Array`, not a `Vector`:
	sig = conv2seriesmatrix2D(sig)

	#Generate matrix of sample numbers:
	M = size(sig,1) #Number of rows (ex: quantizers)
	N = size(sig,2) #Number of samples
	n = 1:N
	sn = ones(Float64, M)*n' #Matrix of "sample number" (Use Float64s for downstream functions)

	#Generate staircase waveform of signal:
	sig = wfrm_stairs(sn, sig)

	push!(plot,
		cons(:wfrm, sig, line=set(style=:solid, color=color, width=2), label=label),
	)
	return plot
end

function plotModTransient(inputSig, outputSig; legend::Bool=true, color=:blue)
	plot = plotModTransient(legend=legend)
	plotModTransient!(plot, outputSig, color=color, label="output")
	plotModTransient!(plot, inputSig, color=:red, label="input")
	return plot
end


"""`lollipop(x, y, color=:blue, lw::Float64=2.0, ybot=0)`

Plot lollipops (o's and sticks).
"""
function plotLollipop(x, y, color=:blue, lw::Float64=2.0, ybot::Float64=0.0, label="", legend=false)
	plot = cons(:plot, linlin, title = "Time-domain", legend=legend,
		labels=set(xaxis="Time", yaxis="Amplitude"),
	)
	simglyph = cons(:a, glyph=set(shape=:o, size=1.5, color=color, fillcolor=color))

	#Make sure we have column vectors of Float64 (Required by InspectDR):
	x = Float64.(x[:]); y = Float64.(abs.(y[:]))
	len = length(x)
	if length(y) != len
		error("x & y array lengths must match.")
	end

	#Generate sticks
	xsticks = [x'; x'; repeat([NaN], len)']
	ysticks = [y'; repeat([ybot], len)'; repeat([NaN], len)']
	xsticks = xsticks[:]; ysticks = ysticks[:] #Need to re-sort as vectors

	#Generate waveforms:
	values = waveform(x, y)
	sticks = waveform(xsticks, ysticks)

	push!(plot,
		cons(:wfrm, sticks, line=set(style=:solid, color=color, width=lw), label=label),
		cons(:wfrm, values, simglyph, line=set(style=:none, color=color, width=lw), label=label),
	)

	return plot
end


#Last line
