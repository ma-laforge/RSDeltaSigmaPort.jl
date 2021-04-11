#RSDeltaSigmaPort: Dataset manipulations
#-------------------------------------------------------------------------------

#Validate x, y pair simple vectors
function validatexy_vec(x, y)
	sz = size(x); sz_y = size(y)
	if sz != sz_y
		throw("x & y sizes do not match")
	elseif length(sz) > 1
		throw("Dimension of data too large (Simple vector expected).")
	end
	return
end

#Validate x, y pair for 1 parameter sweep (2D array)
function validatexy_1param(x, y)
	sz = size(x); sz_y = size(y)
	if sz != sz_y
		throw("x & y sizes do not match")
	elseif length(sz) > 2
		throw("Dimension of data too large (2D arrays expected).")
	end
	return
end

#Last line
