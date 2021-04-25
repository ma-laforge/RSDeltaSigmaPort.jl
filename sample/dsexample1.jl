#Example lowpass/bandpass real/quadrature modulator design
using RSDeltaSigmaPort
using RSDeltaSigmaPort.EasyPlot #set, cons
import Printf: @sprintf
j=im


"""`result = dsexample1(dsm; ampdB::Float64=-3.0, ftest=nothing, LiveDemo=false)`

Example lowpass/bandpass real/quadrature modulator design.

# Output `result`: `Dict` of:
 - :nlev, :NTF, :ABCD, :umax, :amp, :snr, :peak_snr
 - :coeff (a NamedTuple of: a, g, b c)
"""
function dsexample1(dsm::AbstractDSM; ampdB::Float64=-3.0, ftest=nothing, LiveDemo::Bool=false)
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

	LiveDemo && (@info("Realizing NTF..."); flush(stdout); flush(stderr))
#-------------------------------------------------------------------------------
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
	LiveDemo && println("\tdone.")
	plot = documentNTF(dsm, ABCD, sizes=sizes, frespOnly=false)
	saveimage(:png, "dsexample1_NTF.png", plot, AR=2/1, width=900)
	displaygui(plot)

	LiveDemo && (@info("Performing time-domain simulations..."); flush(stdout); flush(stderr))
#-------------------------------------------------------------------------------
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
	LiveDemo && println("\tdone.")
	plot = plotModTransient(u, simresult.v, legend=false)
		ylim = (dsm.M+0.25)
		set(plot, xyaxes=set(ymin=-ylim, ymax=ylim))
	displaygui(plot)

	#Example spectrum
#-------------------------------------------------------------------------------
	plot = plotExampleSpectrum(dsm, NTF, ampdB=ampdB, ftest=ftest, sizes=sizes)
	displaygui(plot)

	#SQNR plot
#-------------------------------------------------------------------------------
	snrinfo = calcSNRInfo(dsm, NTF=NTF)
	plot = plotSNR(snrinfo, dsm)
		#TODO: Apply sizes
	displaygui(plot)

	if isquadrature(dsm) #Example I/Q mismatch
		error("Quadrature section not ported over.")
	end

	umax = nothing
	coeffs = nothing
	if !isquadrature(dsm) #Dynamic range scaling
		LiveDemo && (@info("Performing dynamic range scaling..."); flush(stdout); flush(stderr))
#-------------------------------------------------------------------------------
		ABCD0 = ABCD
		ABCD, umax = scaleABCD(ABCD0, nlev=nlev, f=dsm.f0)
		a,g,b,c = mapABCD(ABCD, dsm.form)
		coeffs = (a=a, g=g, b=b, c=c)
		LiveDemo && println("\tdone.")

		LiveDemo && (@info("Verifying dynamic range scaling..."); flush(stdout); flush(stderr))
#-------------------------------------------------------------------------------
		N = 10^4; N0 = 50
		u = range(0, stop=0.95*umax, length=30)
		test_tone = cos.(2π*dsm.f0 * (0:N-1))
		test_tone[1:N0] = test_tone[1:N0] .* (0.5 .- 0.5*cos.(2π/N0 * (0:N0-1)))

		plot = plotStateMaxima(u, ABCD, test_tone, nlev=nlev)
			set(plot, RSDeltaSigmaPort.linlin, xyaxes=set(ymin=0, ymax=1))
		LiveDemo && println("\tdone.")
		displaygui(plot)
	end

	results = Dict{Symbol, Any}(
		:nlev => nlev,
		:NTF => NTF,
		:ABCD => ABCD,
		:SNR_vs_amp_sim => snrinfo.vs_amp_sim,
		:peak_snr => snrinfo.peak[2],
		:umax => umax,
		:coeffs => coeffs
	)
	return results
end

opt=2
#dsm = RealDSM(order=5, OSR=32, form=:CRFB, Hinf=1.5, opt=2)
dsm = RealDSM(order=5, OSR=32, form=:CRFB, Hinf=1.5, opt=opt)
results = dsexample1(dsm, LiveDemo=true)

:END_OF_EXAMPLE
