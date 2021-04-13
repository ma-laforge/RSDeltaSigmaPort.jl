#RSDeltaSigmaPort: Optimization functions
#-------------------------------------------------------------------------------

"values for optimal zeros when opt==1"
const _VALUES_OPTZEROS_OPT1 = [ #Currently not type stable
	0.0, # n = 1
	sqrt(1/3),
	[sqrt(3/5), 0],
	sqrt.([3/7+sqrt(9/49-3/35), 3/7-sqrt(9/49-3/35)]),
	sqrt.([5/9+sqrt(25/81-5/21), 5/9-sqrt(25/81-5/21), 0]), # n = 5
	[0.23862059, 0.66120988, 0.9324696],
	[0, 0.40584371, 0.74153078, 0.94910785],
	[0.18343709, 0.52553345, 0.79666684, 0.96028993],
	[0, 0.32425101, 0.61337056, 0.83603082, 0.9681602],
	[0.1834370913, 0.5255334458, 0.7966668433, 0.9602899327], # n = 10
	[0, 0.26953955, 0.51909468, 0.73015137, 0.88706238, 0.97822864],
	[0.12523875, 0.36783403, 0.58731921, 0.7699033, 0.90411753, 0.9815607],
	[0, 0.23045331, 0.44849063, 0.64234828, 0.8015776, 0.91759824, 0.98418306],
	[0.10806212, 0.31911586, 0.51525046, 0.68729392, 0.82720185, 0.92843513, 0.98628389], # n = 14
]

"values for optimal zeros when opt==2" #or other?; is it called in other cases?
const _VALUES_OPTZEROS_OPT2 = [
	_VALUES_OPTZEROS_OPT1[1], # n = 1
	0.0,
	_VALUES_OPTZEROS_OPT1[3],
	[0, sqrt(5/7)],
	_VALUES_OPTZEROS_OPT1[5], # n = 5
	sqrt.([0, 7/11+sqrt(56)/33, 7/11-sqrt(56)/33]),
	_VALUES_OPTZEROS_OPT1[7],
	[0, 0.50563161, 0.79017286, 0.95914731],
	_VALUES_OPTZEROS_OPT1[9],
	[0, 0.41572267, 0.67208682, 0.86238894, 0.97342121], # n = 10
	_VALUES_OPTZEROS_OPT1[11],
	[0, 0.35222363, 0.58006251, 0.76647993, 0.90281326, 0.98132047],
	_VALUES_OPTZEROS_OPT1[13],
	[0, 0.30524384, 0.50836649, 0.6836066, 0.82537239, 0.92772336, 0.98615167], # n = 14
]
@assert(length(_VALUES_OPTZEROS_OPT1) == length(_VALUES_OPTZEROS_OPT2))

"""`optZeros = ds_optzeros( n, opt=1 )`
A helper function for the synthesizeNTF function of the Delta-Sigma Toolbox.
Returns the zeros which minimize the in-band noise power of 
a delta-sigma modulator's NTF.
"""
function ds_optzeros(n, opt=1)
	norig = n; n = array_round(n) #VERIFYME
	local optZeros

	if opt==0
		optZeros = zeros(ceil(Int, norig/2))
	elseif n > length(_VALUES_OPTZEROS_OPT1)
		throw(ErrorException("Optimized zeros for n>14 are not available."))
	elseif opt==1
		optZeros = _VALUES_OPTZEROS_OPT1[n]
	else
		optZeros = _VALUES_OPTZEROS_OPT2[n]
	end

	# Sort the zeros and replicate them.
	z = sort(optZeros)
	optZeros = zeros(n)
	m=1

	if mod(n,2)==1
		optZeros[1] = z[1]
		z = z[2:length(z)]
		m=m+1
	end
	for i=1:length(z)
		optZeros[m]   =  z[i]
		optZeros[m+1] = -z[i]
		m = m+2
	end
	return optZeros
end

#Last line
