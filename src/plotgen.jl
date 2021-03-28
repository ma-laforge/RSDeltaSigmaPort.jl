#RSDeltaSigmaPort: Base plot generation tools
#-------------------------------------------------------------------------------

#==Misc. draw functions
===============================================================================#
function plot_unitcircle!(plot; color=:black, width=2, nseg::Int=360)
	Θ = range(0, stop=2π, length=nseg+1)
	cplxcircle = exp.((2π*j)*Θ)
	circle = DataF1(real(cplxcircle), imag(cplxcircle))
	push!(plot,
		cons(:wfrm, circle, line=set(style=:solid, color=color, width=width))
	)
	return plot
end


#==Pole-zero plots
===============================================================================#
#Create empty pole-zero plot
function plotPZ()
	plot = cons(:plot, title = "Pole-Zero Plot", legend=false,
		xyaxes=set(xscale=:lin, yscale=:lin, xmin=-1.1, xmax=1.1, ymin=-1.1, ymax=1.1),
	)
	plot_unitcircle!(plot)
	return plot
end

"""`function plotPZ!(H,color=:black,markersize=5,list=false)`
Plot the poles and zeros of a transfer function.
If list=true, a list of the poles and zeros is superimposed on the plot.
"""
function plotPZ!(plot, H::ZPKData; color=:blue, markersize=2, list::Bool=false)
	isVector = isa(color, Vector)
	color_pole = isVector ? color[1] : color
	color_zero = isVector ? color[end] : color
	attrib_poles = cons(:a,
		glyph=set(shape=:x, color=color_pole, size=markersize), line=set(style=:none, width=2)
	)
	attrib_zeros = cons(:a,
		glyph=set(shape=:o, color=color_zero, size=markersize), line=set(style=:none, width=2)
	)
	
	z, p, k = _zpkdata(H)

	# Plot x and o for poles and zeros, respectively:
	push!(plot, cons(:wfrm, DataF1(real(p), imag(p)), attrib_poles))
	if !isempty(z) #Add zeros
		push!(plot, cons(:wfrm, DataF1(real(z), imag(z)), attrib_zeros))
	end

#=
	if list
		# List the poles and zeros
		p = cplxpair(p);
		y = 0.05*(ceil(length(p)/2)+1);
		str = 'Poles:               ';
		text( 0, y, str, 'Hor', 'Right', 'Ver', 'Mid'); y = y - 0.1;
		for i = 1:2:length(p);
			if abs(imag(p(i))) < 1e-6
				str = sprintf('%+.4f      ', real(p(i)) );
			else
				str = sprintf( '%+.4f+/-j%.4f  ', real(p(i)), abs(imag(p(i))) );
			end
			text( 0, y, str, 'Hor', 'Right', 'Ver', 'Mid'); y = y - 0.1;
		end
		if !isempty(z)
			z = z( !isnan(z) && !isinf(z) );
			z = cplxpair(z);
			y = 0.05*(ceil(length(z)/2)+1);
			str = '        Zeros:';
			text( 0, y, str, 'Hor', 'Left', 'Ver', 'Mid'); y = y - 0.1;
			for i = 1:2:length(z);
				if abs(imag(z(i))) < 1e-6
					str = sprintf('%+.4f      ', real(z(i)) );
				else
					str = sprintf( '  %+.4f+/-j%.4f', real(z(i)), abs(imag(z(i))) );
				end
				text( 0, y, str, 'Hor', 'Left', 'Ver', 'Mid');
				y = y - 0.1;
			end
		end
	end
=#
	return plot
end

plotPZ(H::ZPKData, args...; kwargs...) = plotPZ!(plotPZ(), H::ZPKData, args...; kwargs...)


#==NTF plots
===============================================================================#
function plotNTF()
	dBmin = -100 #Limit dB window to finite value
	𝑓min_log = 10^-2 #Limit min 𝑓 for log(f) plots
	nogrid = set(vmajor = false, vminor = false, hmajor=false, hminor=false)
	vgridonly = set(vmajor = true, vminor = true, hmajor=false, hminor=false)

	plot_PZ = plotPZ()
	plot_PZ.title = "NTF Poles and Zeros"

	plot_NTFlin = cons(:plot, legend=false,
		title = "NTF Magnitude Response",
		xyaxes=set(xscale=:lin, yscale=:lin, ymin=dBmin),
		labels=set(xaxis="Normalized frequency (1→f_s)", yaxis="dB"),
	)

	plot_NTFlog = cons(:plot, legend=false, grid = vgridonly,
#		title = "NTF Magnitude Response",
		xyaxes=set(xscale=:log, yscale=:lin, xmin=𝑓min_log, ymin=dBmin),
		labels=set(xaxis="Normalized frequency (1→f_B)", yaxis="dB"),
	)

	pcoll = push!(cons(:plot_collection, title="Sample Plot"),
		plot_PZ, plot_NTFlin, plot_NTFlog
	)

	pcoll.bblist = [
		BoundingBox(0, 0.5, 0, 1), #Pole-zero
		BoundingBox(0.5, 1, 0, 0.5), #NTF linf
		BoundingBox(0.5, 1, 0.5, 1), #NTF logf
	]

	return pcoll
end

function plotNTF!(pcoll, NTF::ZPKData, OSR::Int; color=:blue)
	dfltline = cons(:a, line=set(style=:solid, color=color, width=3))
	markerline = cons(:a, line=set(style=:dash, width=2.5))
	plot_PZ, plot_NTFlin, plot_NTFlog = pcoll.plotlist

	plotPZ!(plot_PZ, NTF, color=color)

	𝑓 = vcat(range(0, stop=0.75/OSR, length=100), range(0.75/OSR, stop=0.5, length=100))
	z = exp.(2j*π*𝑓)
	magNTF = dbv(evalTF(NTF, z))
	push!(plot_NTFlin, 
		cons(:wfrm, DataF1(𝑓, magNTF), dfltline, label="|NTF|")
	)

	𝑓start = 0.01
	𝑓 = collect(range(𝑓start, stop=1.2, length=200)/(2*OSR))
	z = exp.(2j*π*𝑓)
	𝑓norm = 𝑓*(2*OSR)
	magNTF = dbv(evalTF(NTF, z))
	σ_NTF = dbv(rmsGain(NTF, 0, 0.5/OSR))
	rmsgain_str = @sprintf("RMS Gain = %5.0fdB", σ_NTF)
	push!(plot_NTFlog, 
		cons(:wfrm, DataF1(𝑓norm, magNTF), dfltline, label="|NTF|"),
		cons(:atext, rmsgain_str, y=σ_NTF, offset=set(y=3), reloffset=set(x=0.5), align=:bc),
		cons(:hmarker, σ_NTF, markerline, strip=2),
	)

	return pcoll
end

plotNTF(NTF::ZPKData, args...; kwargs...) = plotNTF!(plotNTF(), NTF::ZPKData, args...; kwargs...)


#Last line
