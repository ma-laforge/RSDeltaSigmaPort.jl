#RSDeltaSigmaPort: Generate plots of the Signal-to-Noise Ratio
#-------------------------------------------------------------------------------

#==SNR plot
===============================================================================#
#TODO: Rename `plotSQNR`
function plotSNR(; legend::Bool=true,
		title::String="SNR: Theory vs Simulation"
	)

	plot = cons(:plot, linlin, title=title, legend=legend,
		#xyaxes=set(xmin=-120, xmax=0, ymin=0, ymax=120),
		xyaxes=set(ymax=120), #Clip infinite gains only
		labels=set(xaxis="Input Level (dBFS)", yaxis="SQNR (dB)"),
	)

	return plot
end

function plotSNR!(plot, sig, NTF, OSR::Int; f0::Float64=0.0, nlev::Int=2,
		id::String="simulation", color=:green
	)
	sig = conv2seriesmatrix2D(sig)
	M, N = size(sig)

	snr, amp = simulateSNR(NTF, OSR, f0=f0, nlev=nlev)
	pk_snr, pk_amp = peakSNR(snr, amp) #Single values

	snr_pred, amp_pred = nothing, nothing
	if nlev == 2 #predictSNR() only supports 2 levels
		snr_pred, amp_pred = predictSNR(NTF, OSR, f0=f0)
		amp_pred = conv2seriesmatrix2D(amp_pred*1.0) #Ensure amplitudes are Float64
	end

	#Convert to waveforms:
	snr = waveform(amp, snr)

	simglyph = cons(:a, glyph=set(shape=:o, size=1, color=color, fillcolor=color))

	snr_str = @sprintf("peak SNR = %4.1fdB\n@ A = %4.1f dBFS (OSR = %d)",
		pk_snr, pk_amp, OSR #orig. coords: (-25, 85)
	)
	push!(plot, 
		cons(:wfrm, snr, simglyph, line=set(style=:dashdot, color=color, width=2), label="simulation"),
		cons(:atext, snr_str, x=-25, y=85, align=:cr),
	)
	if !isnothing(snr_pred)
		snr_pred = waveform(amp_pred, snr_pred) #Convert to waveform
		push!(plot,
			cons(:wfrm, snr_pred, line=set(style=:solid, color=:red, width=2), label="theory"),
		)
	end
	return plot
end

function plotSNR(sig, NTF, OSR::Int; f0::Float64=0.0, nlev::Int=2, legend::Bool=true,
		title::String="SQNR vs. Input Level", id::String="simulation", color=:green
	)
	plot = plotSNR(legend=legend, title=title)
	plotSNR!(plot, sig, NTF, OSR, f0=f0, nlev=nlev, id=id, color=color)
	return plot
end


#Last line
