#RSDeltaSigma: Tools used to aid with simulations
#-------------------------------------------------------------------------------


#==selectElement()
===============================================================================#
"""`sv = selectElement(v, sy, dw=[], tri=false)`

Select elements of a multi-element DAC to minimize the selection error, 
subject to the constraint that the nominal DAC output is v, i.e. v = sv'*dw.
Assume that the preferred usage order is that given by dw.

# Input
 - `v`:   DSM output in `-M:2:M`
 - `sy`:  Mx1 desired usage vector
 - `dw`:  Mx1 vector of element weights
 - `tri`: Flag indicating the elements are tri-level, i.e. `sv=0` is allowed.
   NOT SUPPORTED YET.

# Output
 - `sv`: Mx1 selection vector. `+/-1` if `tri==false`, `{0,+-1}` if `tri==true`.
"""
function selectElement(v, sy::Vector, dw::Vector=[], tri=false) %#ok<INUSD>
	M = length(sy)
	sdw = M
	dw_is_one = (isempty(dw) || all(dw .== 1))
	if !dw_is_one
		sdw = sum(dw)
	end

	if abs(v) > sdw
		throw("|v| is too large.")
	elseif v == sdw
		return ones(M,1)
	elseif v == -sdw
		return -ones(M,1)
	end

	sv = -ones(M,1)
	possibilities = sortperm(-sy) #Determine the element priority

	#Speed up selection in the usual case where all element weights are one.
	#Suggested 1997-03-05 by J. A. Cherry.
	if dw_is_one
		sv[possibilities[1:array_round((M+v)/2)]] .= 1
		return sv
	end

	#Go through sv possibilities one by one, until one which meets the
	#v = sv'*dw constraint is found.
	i = 1 #Selection level
	pointer = ones(1,M) #Array of pointers to selected elements
	selected = zeros(1,M) #Selected elements
	dw0 = [0; dw]

	while true
		while pointer[i]>M
			#backtrack
			i = i-1
			if i==0
				break #failure!
			end
			pointer[i] = pointer[i] + 1
			selected[i+1:end] = 0
		end
		selected[i] = possibilities[pointer[i]]
		dv = 2*sum(dw0[selected+1]) - sdw
		if dv==v
			break #success!
		elseif dv<v
			#Proceed to the next level of selection
			pointer[i+1] = pointer[i] + 1
			i = i+1
		else
			#Try the next element at the current level.
			pointer[i] = pointer[i] + 1
		end
	end
	selected = selected[findall(selected .!= 0)]
	sv[selected] = 1
	return sv
end

#Last line
