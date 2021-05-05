#RSDeltaSigmaPort: Calculate frequency spectrum information
#-------------------------------------------------------------------------------


#==calcSpecInfo
===============================================================================#
"""
Calculate spectrum information

# Inputs
 - sig: Signal to analyze
 - NTF: Noise transfer function to compute theoretical noise spectrum
 - fband: Specifies (𝑓min, 𝑓max) band limits of the modulator (normalized 𝑓)
 - ftest: test frequency (normalized 𝑓)
 - M: Number of quantizer steps (nlev = M+1)
 - quadrature:
"""
function calcSpecInfo(sig, NTF, fband, ftest::Float64; M::Int=1, quadrature::Bool=false)
	sig = conv2seriesmatrix2D(sig)
	R, N = size(sig)
	if R!=1
		throw("FIXME: replicate hann window rows for R>0")
	end

	iband = round(Int, minimum(fband)*N):round(Int, maximum(fband)*N)
	itest = round(Int, ftest*N)

	Δ = 2
	NBW = 1.5/N
	spec0 = fft(applywnd(sig, ds_hann(N)))/(M*N/4)
	local 𝑓, specdB, SNR, ePSD
	if !quadrature
		𝑓pts = div(N,2)+1 #Only need half the spectrum for display purposes
		𝑓 = range(0, stop=0.5, length=𝑓pts)
		specdB = dbv.(spec0[1:𝑓pts])
#		@show iband .+ 1, itest-minimum(iband)
		SNR = calculateSNR(spec0[iband .+ 1], itest-minimum(iband))

		#Compute Expected power spetral density:
#		Sqq = 4 * evalTF(NTF ,exp.(2π*j*𝑓)).^2 / (3*M^2)
		Snn = abs.(evalTF(NTF, exp.(2π*j*𝑓))).^2 * 2*2/12 * (Δ/M)^2
		ePSD = dbp.(Snn*NBW)
	else
		throw("Not implemented")
	end

	return (freq=𝑓, spec0=spec0, specdB=specdB, ePSD=ePSD,
		SNR=SNR, NBW=NBW, NTF=NTF, M=M, iband=iband, itest=itest
	)
end

function calcSpecInfo(sig, dsm::RealDSM, fband=nothing, ftest=nothing)
	if isnothing(fband); fband = default_fband(dsm); end
	if isnothing(ftest); ftest = default_ftest(dsm); end
	NTF = synthesizeNTF(dsm)
	return calcSpecInfo(sig, NTF, fband, ftest, M=dsm.M, false)
end

#Last line
