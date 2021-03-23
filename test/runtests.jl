using Test, RSDeltaSigmaPort


function printheader(title)
	println("\n", title, "\n", repeat("=", 80))
end

function show_testset_section()
	println()
	@info "SECTION: " * Test.get_testset().description * "\n" * repeat("=", 80)
end

function show_testset_description()
	@info "Testing: " * Test.get_testset().description
end

testfiles = ["mlfunc.jl"]

for testfile in testfiles
	include(testfile)
end

:Test_Complete
