#RSDeltaSigma: Other math functions
#-------------------------------------------------------------------------------

"""`A  = mapQtoR(Z)`

Convert a quadrature matrix into its real equivalent.

Each element in `Z` is represented by a 2x2 matrix in `A`:
 - `z -> [x -y; y x]`
"""
function mapQtoR(Z)
	A = zeros(2*size(Z))
	A[1:2:end,1:2:end] = real(Z)
	A[2:2:end,2:2:end] = real(Z)
	A[1:2:end,2:2:end] = -imag(Z)
	A[2:2:end,1:2:end] = imag(Z)
	return A
end

#Last line
