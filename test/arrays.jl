@testset "Array-related tests" begin show_testset_section() #Scope for test data

@testset "pad*()" begin show_testset_description()
	v = [4,5,6,10]
	@test padt(v, 7) == [0,0,0, 4,5,6,10]
	@test padl(v, 7) == [0,0,0, 4,5,6,10]
	@test padb(v, 7) == [4,5,6,10, 0,0,0]
	@test padr(v, 7) == [4,5,6,10, 0,0,0]
	v = [4 5 6; 8 10 11]
	@test padt(v, 5) == [0 0 0; 0 0 0; 0 0 0; 4 5 6; 8 10 11]
	@test padb(v, 5) == [4 5 6; 8 10 11; 0 0 0; 0 0 0; 0 0 0]
	@test padl(v, 7) == [0 0 0 0 4 5 6; 0 0 0 0 8 10 11]
	@test padr(v, 7) == [4 5 6 0 0 0 0; 8 10 11 0 0 0 0]
end

end
