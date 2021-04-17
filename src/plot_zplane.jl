#RSDeltaSigmaPort: Generate plots of the complex Z-plane
#-------------------------------------------------------------------------------


#==Waveform builders
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
		labels=set(xaxis="Real", yaxis="Imag"),
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


#Last line
