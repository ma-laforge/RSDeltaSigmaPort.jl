#Example lowpass/bandpass real/quadrature modulator design
using RSDeltaSigmaPort
using RSDeltaSigmaPort: AbstractDSM, isquadrature
using RSDeltaSigmaPort.EasyPlot #set, cons
import Printf: @sprintf
j=im



"""`result = dsexample1(dsm, Atest::Float64=-3.0, Ftest=nothing, LiveDemo=false)`

Example lowpass/bandpass real/quadrature modulator design.

# Output `result`: `Dict` of:
 - :nlev, :ntf, :ABCD, :umax, :amp, :snr, :peak_snr
 - :coeff (a NamedTuple of: a, g, b c)
"""
function dsexample1(dsm::AbstractDSM, Atest::Float64=-3.0, Ftest=nothing, LiveDemo::Bool=false)
	#Computed defaults
	if nothing == Ftest
		if 0==dsm.f0 || isquadrature(dsm)
			Ftest = 0.15/dsm.osr
		else
			Ftest = dsm.f0 + 0.08/dsm.osr
		end
	end

	#Derived parameters
	nlev = dsm.M + 1
	typestr = (0==dsm.f0) ? "Lowpass" : "Bandpass"
	if isquadrature(dsm)
		typestr = "Quadrature" * typestr
	end

	println("\t", ds_orderString(dsm.order,true), "-Order $typestr Example... ")
	#TODO: support plot attributes?:
	sizes = Dict{Symbol, Any}(
		:lw => 1, #LineWidth
		:ms => 5, #MarkerSize
		:fs => 12, #FontSize
		:fw => "normal", #FontWeight
	)
	if LiveDemo #Overwrite settings:
		push!(sizes,
			:lw => 2,
			:ms => 6,
			:fw => "bold",
		)
	end

	#NTF synthesis and realization
	#fprintf(1,'Doing NTF synthesis and realization... ');
	local ntf, ABCD
	if !isquadrature(dsm)
		ntf = synthesizeNTF(dsm)
		a,g,b,c = realizeNTF(ntf, dsm.form)
		z0 = exp(2π*j*dsm.f0)
		b = [abs(b[1]+b[2]*(1-z0)) zeros(1,length(b)-1)] #Use a single feed-in for the input
		ABCD = stuffABCD(a,g,b,c,dsm.form)
	else
		ntf = synthesizeQNTF(dsm)
		ABCD = realizeQNTF(ntf, dsm.form, 1)
	end
	#fprintf(1,'Done.\n');
	plot = documentNTF(ABCD, dsm.osr, dsm.f0, isquadrature(dsm), sizes, 1)
#	displaygui(plot)

	#Time-domain simulations
	#fprintf(1,'Doing time-domain simulations... ');

	#Time-domain plot
	N = 100
	t = 0:N-1
	local u, v
	if !isquadrature(dsm)
		u = undbv(Atest)*dsm.M*sin.( 2π*Ftest*t )
		(v,) = simulateDSM(u, ABCD, nlev = nlev)
	else
		u = undbv(Atest)*dsm.M*exp.( 2π*j*Ftest*t )
		(v,) = simulateQDSM(u, ABCD, nlev=nlev)
	end
	plot = plotModTransient(u, v, legend=false)
		#CHECK FOR COLORS
	ypk = (dsm.M+0.25)
	set(plot, xyaxes=set(ymin=-ypk, ymax=ypk))
	displaygui(plot)

	#Example spectrum
	plot = plotExampleSpectrum(dsm, sizes=sizes, Atest=Atest, Ftest=Ftest)
	displaygui(plot)

	#SQNR plot
	plot = plotSNR(v, ntf, dsm.osr, f0=dsm.f0)
		#TODO: Apply sizes
	displaygui(plot)

	results = Dict{Symbol, Any}(
		:nlev => nlev,
		:ntf => ntf,
		:ABCD => ABCD,
#		:amp => amp,
#		:snr => snr,
#		:umax => umax,
#		:peak_snr => peak_snr,
	)

	if !LiveDemo && !isquadrature(dsm)
		push!(results, :coeffs=>(a=a, g=g, b=b, c=c))
	end

	return results
end

dsm = RealDSM()
results = dsexample1(dsm)


#Last line
