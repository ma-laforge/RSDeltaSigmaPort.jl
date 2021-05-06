#RSDeltaSigma: exampleHBF() function
#-------------------------------------------------------------------------------

function HBFex1()
	q1 = [
		 0 -4 -5;
		 1  1  1;
		 0 -4 -5;
		-1 -1 -1;
		-1 -3 -5;
		 1  1  1;
		-3 -5  0;
		-1 -1  0;
	]
	q2 = [
		-1 -3 -7;
		 1  1  1;
		-2 -5 -7;
		-1  1  1;
		-3 -8  0;
		 1 -1  0;
		-4 -6 -8;
		-1 -1 -1;
		-4 -8  0;
		 1 -1  0;
		-5 -7 -8;
		-1 -1 -1;
		-5  0  0;
		 1  0  0;
		-6 -7  0;
		-1 -1  0;
		-7 -8  0;
		-1 -1  0;
		-6  0  0;
		 1  0  0;
	]
	return (q1,q2)
end

"""`HBFex2()`

Coefficients from toolbox v2.0 documentation.
"""
function HBFex2()
	q1 = [
		 0 -4 -7;
		 1 -1  1;
		-1 -3 -6;
		-1 -1 -1;
		-2 -4 -7;
		 1 -1  1;
	]
	q2 = [
		-1 -3 -8;
		 1  1 -1;
		-2 -4 -9;
		-1  1 -1;
		-3 -5 -9;
		 1 -1  1;
		-4 -7 -8;
		-1  1  1;
		-5 -8 -11;
		 1 -1 -1;
		-6 -9 -11;
		-1  1 -1;
	]
	return (q1,q2)
end

"""`Saramaki88()`

Coefficients from "Efficient VLSI-Realizable Decimators for Sigma-Delta
Analog-to-Digital Converters," T. Saramaki and H. Tenhunen, ISCAS 1988,
pp 1525-2528.
"""
function Saramaki88()
	q1 = [
		 0 -2;
		 1 -1;
		-2 -9;
		-1  1;
	]
	q2 = [
		-1 -3  0;
		 1  1  0;
		-3 -5 -6;
		-1 -1 -1;
		-4 -6  0;
		 1  1  0;
	]
	return (q1,q2)
end

"""`Saramaki90()`

Coefficients from "Multiplier-Free Decimation Algorithms for Superresolution
Oversampled Converters," T. Saramaki, T. Karema, T. Ritoniemi, and H. Tenhunen
ISCAS 1990, vol. 4. pp 1525-2528.
"""
function Saramaki90()
	q1 = [
		 0 -4 -5;
		 1  1  1;
		 0 -4 -5;
		-1 -1 -1;
		-1 -3 -5;
		 1  1  1;
		-3 -5 -20;
		-1 -1  1;
	]
	q2 = [
		-1 -3 -7;
		 1  1  1;
		-2 -5 -7;
		-1  1  1;
		-3 -8  0;
		 1 -1  0;
		-4 -6 -8;
		-1 -1 -1;
		-4 -8  0;
		 1 -1  0;
		-5 -7 -8;
		-1 -1 -1;
		-5  0  0;
		 1  0  0;
		-6 -7  0;
		-1 -1  0;
		-6  0  0;
		 1  0  0;
		-7 -8  0;
		-1 -1  0;
		-6  0  0;
		 1  0  0;
	]
	return (q1,q2)
end

function exampleHBF(n::Int = 0)
	fnlist = [HBFex1, HBFex2, Saramaki88, Saramaki90]
	if !in(n, keys(fnlist))
		throw("exampleHBF: No code for example $n.")
	end
	fn = fnlist[n]
	(q1, q2) = fn()
	F1 = convertForm(bunquantize(q1),q1)
	F2 = convertForm(bunquantize(q2),q2)
	return (F1, F2)
end

#Last line
