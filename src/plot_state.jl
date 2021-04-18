#RSDeltaSigmaPort: Generate plots of internal state
#-------------------------------------------------------------------------------


#==Plot state maxima
===============================================================================#
function plotStateMaxima(u, ABCD; nlev::Int=2, N::Int=10^5)
	title = "Simulated State Maxima"; color = :blue
	T = ones(1,N)
	nu, nq, order = get_nu_nq_order(T, ABCD, nlev)
	u = collect(u)
	maxima = zeros(order,length(u))

	for i = 1:length(u)
		ui = u[i]

		simresult = simulateDSM(ui*T, ABCD, trackmax=true)
		maxima[:,i] = simresult.xmax[:]
		if any(simresult.xmax .> 1e2)
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
	for i in 1:order
		id = "state $i"
		_maxima = waveform(u, maxima[i,:])
		push!(plot,
			cons(:wfrm, _maxima, line=set(style=:solid, width=2), label=id),
		)
	end
	return plot
end

#Last line
