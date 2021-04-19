#RSDeltaSigmaPort: Generate plots of the Signal-to-Noise Ratio
#-------------------------------------------------------------------------------

#==SNR plot
===============================================================================#
#TODO: Rename `plotSQNR`?
function plotSNR()
	plot = cons(:plot, linlin, title="SQNR vs. Input Level", legend=true,
		#xyaxes=set(xmin=-120, xmax=0, ymin=0, ymax=120),
		xyaxes=set(ymax=120), #Clip infinite gains only
		labels=set(xaxis="Input Level (dBFS)", yaxis="SQNR (dB)"),
	)
	return plot
end

function plotSNR!(plot, snrinfo, dsm::RealDSM; id::String="simulation", color=:green)
	simglyph = cons(:a, glyph=set(shape=:o, size=1, color=color, fillcolor=color))

	#Access data:
	(amp, SNR) = snrinfo.vs_amp_sim
	(pk_amp, pk_SNR) = snrinfo.peak

	#Convert to waveforms:
	SNR = waveform(amp, SNR)

	infostr = @sprintf("peak SNR = %4.1fdB\n@ A = %4.1f dBFS (OSR = %d)",
		pk_SNR, pk_amp, dsm.OSR #orig. coords: (-25, 85)
	)
	push!(plot, 
		cons(:wfrm, SNR, simglyph, line=set(style=:dashdot, color=color, width=2), label="simulation"),
		cons(:atext, infostr, x=-25, y=85, align=:cr),
	)
	if !isnothing(snrinfo.vs_amp_predicted)
		(amp_pred, SNR_pred) = snrinfo.vs_amp_predicted
		SNR_pred = waveform(amp_pred, SNR_pred) #Convert to waveform
		push!(plot,
			cons(:wfrm, SNR_pred, line=set(style=:solid, color=:red, width=2), label="theory"),
		)
	end
	return plot
end

function plotSNR(snrinfo, dsm::RealDSM; id::String="simulation", color=:green)
	plot = plotSNR()
	plotSNR!(plot, snrinfo, dsm, id=id, color=color)
	return plot
end


#Last line
