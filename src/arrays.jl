#RSDeltaSigma: Array manipulations
#-------------------------------------------------------------------------------


#==Converters to make input data uniform
===============================================================================#
"""`conv2seriesmatrix2D(x)`

Ensure we have a 2D matrix representing series data (Vector->2×N Array)
"""
conv2seriesmatrix2D(x::T) where T =
	throw(ErrorException("Cannot convert $T to a \"series data matrix\""))
conv2seriesmatrix2D(x::AbstractVector) = collect(x')
conv2seriesmatrix2D(x::Array{T, 2}) where T<:Number = x #Assume format is ok.


#==Array/Vector padding
===============================================================================#
"""`y = padl(x, n, val)`

Pad a matrix x on the left to length n with value val(0)
The empty matrix is assumed to be have 1 empty row
"""
padl(x, n::Int, val=0) =
	[ val*ones(eltype(x), max(1,size(x,1)), n-size(x,2)) x ]
padl(x::Vector, n::Int, val=0) = vcat(val*ones(eltype(x), n-length(x)), x)

"""`y = padr(x, n, val)`

Pad a matrix x on the right to length n with value val(0)
The empty matrix is assumed to be have 1 empty row
"""
padr(x, n::Int, val=0) =
	[ x val*ones(eltype(x), max(1,size(x,1)), n-size(x,2)) ]
padr(x::Vector, n::Int, val=0) = vcat(x, val*ones(eltype(x), n-length(x)))

"""`y = padt(x, n, val)`

Pad a matrix x on the top to length n with value val(0)
The empty matrix is assumed to be have 1 empty column
"""
padt(x, n::Int, val=0) =
	[ val*ones(eltype(x), n-size(x,1), max(1,size(x,2))); x ]
padt(x::Vector, n::Int, val=0) = padl(x,n,val)


"""`y = padb(x, n, val)`

Pad a matrix x on the bottom to length n with value val(0)
The empty matrix is assumed to be have 1 empty column
"""
padb(x, n::Int, val=0) =
	[ x ; val*ones(eltype(x), n-size(x,1), max(1,size(x,2))) ]
padb(x::Vector, n::Int, val=0) = padr(x,n,val)


#==Quadrature <--> real
===============================================================================#
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
