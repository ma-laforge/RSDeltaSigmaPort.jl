#Example lowpass/bandpass real/quadrature modulator design
using RSDeltaSigmaPort
using RSDeltaSigmaPort: AbstractDSM, isquadrature
using RSDeltaSigmaPort.EasyPlot #set, cons
import Printf: @sprintf
j=im



"""`result = dsexample1(dsm, ampdB::Float64=-3.0, ftest=nothing, LiveDemo=false)`

Example lowpass/bandpass real/quadrature modulator design.

# Output `result`: `Dict` of:
 - :nlev, :NTF, :ABCD, :umax, :amp, :snr, :peak_snr
 - :coeff (a NamedTuple of: a, g, b c)
"""
function dsexample1(dsm::AbstractDSM, ampdB::Float64=-3.0, ftest=nothing, LiveDemo::Bool=false)
	#Computed defaults
	if isnothing(ftest)
		ftest = default_ftest(dsm)
	end

	#Derived parameters
	nlev = dsm.M + 1
	dsmtypestr = str_modulatortype(dsm)
	println("\t", dsmtypestr, " Example... ")
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
	local NTF, ABCD
	if !isquadrature(dsm)
		NTF = synthesizeNTF(dsm)
		a,g,b,c = realizeNTF(NTF, dsm.form)
		z0 = exp(2π*j*dsm.f0)
		b = [abs(b[1]+b[2]*(1-z0)) zeros(1,length(b)-1)] #Use a single feed-in for the input
		ABCD = stuffABCD(a,g,b,c,dsm.form)
	else
		NTF = synthesizeQNTF(dsm)
		ABCD = realizeQNTF(NTF, dsm.form, 1)
	end
	#fprintf(1,'Done.\n');
	plot = documentNTF(ABCD, dsm.OSR, dsm.f0, isquadrature(dsm), sizes, 1)
#	displaygui(plot)

	#Time-domain simulations
	#fprintf(1,'Doing time-domain simulations... ');

	#Time-domain plot
	N = 100
	t = 0:N-1
	#Do not use genTestTone() (N too small; ftest gets rounded to 0Hz).
	#(u, iftest) = genTestTone(dsm, ampdB, ftest, N=100)
	local simresult
	if !isquadrature(dsm)
		u = undbv(ampdB)*dsm.M*sin.( 2π*ftest*t )
		simresult = simulateDSM(u, ABCD, nlev=nlev)
	else
		u = undbv(ampdB)*dsm.M*exp.( 2π*j*ftest*t );
		simresult = simulateQDSM(u, ABCD, nlev=nlev)
	end
	plot = plotModTransient(u, simresult.v, legend=false)
		#CHECK FOR COLORS
	ylim = (dsm.M+0.25)
	set(plot, xyaxes=set(ymin=-ylim, ymax=ylim))
	displaygui(plot)

	#Example spectrum
	plot = plotExampleSpectrum(dsm, NTF, ampdB=ampdB, ftest=ftest, sizes=sizes)
	displaygui(plot)

	#SQNR plot
	plot = plotSNR(simresult.v, NTF, dsm.OSR, f0=dsm.f0)
		#TODO: Apply sizes
	displaygui(plot)

	results = Dict{Symbol, Any}(
		:nlev => nlev,
		:NTF => NTF,
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

opt=2
#dsm = RealDSM(order=5, OSR=32, form=:CRFB, Hinf=1.5, opt=2)
dsm = RealDSM(order=5, OSR=32, form=:CRFB, Hinf=1.5, opt=opt)
results = dsexample1(dsm)


#Last line
