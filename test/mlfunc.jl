@testset "Ml function tests" begin show_testset_section() #Scope for test data

p = exp.((2*pi*im/7) .* collect(0:6))

#=
                  1.0 + 0.0im
   0.6234898018587336 + 0.7818314824680298im
 -0.22252093395631434 + 0.9749279121818236im
   -0.900968867902419 + 0.43388373911755823im
  -0.9009688679024191 - 0.433883739117558im
  -0.2225209339563146 - 0.9749279121818236im
   0.6234898018587334 - 0.7818314824680299im
=#

@testset "cplxpair()" begin show_testset_description()
	rtol = 1e-12

	pA = exp.((2*pi*im/7) .* collect(0:6))
	pA = RSDeltaSigmaPort.cplxpair(pA)

	p = -0.900968867902419 + 0.43388373911755823im
		@test pA[1] ≈ conj(p) rtol = rtol
		@test pA[2] ≈ p       rtol = rtol
	p = -0.22252093395631434 + 0.9749279121818236im
		@test pA[3] ≈ conj(p) rtol = rtol
		@test pA[4] ≈ p       rtol = rtol
	p = 0.6234898018587336 + 0.7818314824680298im
		@test pA[5] ≈ conj(p) rtol = rtol
		@test pA[6] ≈ p       rtol = rtol
	p = 1.0 + 0.0im
		@test pA[7] ≈ p rtol = rtol
end
@testset "poly()" begin show_testset_description()
	@test RSDeltaSigmaPort.poly([3,4,5]) == [1, -12, 47, -60]
end

end
