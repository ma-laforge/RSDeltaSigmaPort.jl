#RSDeltaSigmaPort: Generate plots of internal state
#-------------------------------------------------------------------------------


#==Plot state maxima
===============================================================================#
"""`plotStateMaxima(u, ABCD, test_sig=nothing; nlev::Int=2, N::Int=10^5)`

# Inputs
 - `N`: Number of points used to generate a test signal
   (only used if `test_sig=nothing`).
"""
function plotStateMaxima(u, ABCD, test_sig=nothing; nlev::Int=2, N::Int=10^5)
	title = "Simulated State Maxima"; color = :blue
	u = collect(u)
	if isnothing(test_sig)
		test_sig = ones(1,N)
	end
	nu, nq, order = get_nu_nq_order(test_sig, ABCD, nlev)
	maxima = zeros(order,length(u))

	for i = 1:length(u)
		ui = u[i]

		simresult = simulateDSM(ui*test_sig, ABCD, nlev=nlev, trackmax=true)
		maxima[:,i] = simresult.xmax[:]
		if any(simresult.xmax .> 1e2)
			@warn("plotStateMaxima(): umax was too high.")
			umax = ui
			u = u[1:i]
			maxima = maxima[:,1:i]
			break
		end
	end

	plot = cons(:plot, linlog, title = title, legend=true,
#		xyaxes=set(xmin=0, xmax=0.6, ymin=1e-4, ymax=10),
		xyaxes=set(ymin=1e-4, ymax=10),
		labels=set(xaxis="DC input", yaxis="Peak value"),
	)
	colors = range(HSV(0,1,1), stop=HSV(-360,1,1), length=order)
	for i in 1:order
		id = "state $i"; c = colors[i]
		simglyph = cons(:a, glyph=set(shape=:o, size=1, color=c, fillcolor=c))

		_maxima = waveform(u, maxima[i,:])
		push!(plot,
		cons(:wfrm, _maxima, simglyph, line=set(style=:dash, color=c, width=2), label=id),
		)
	end
	return plot
end

#Last line
