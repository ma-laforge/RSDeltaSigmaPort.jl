#RSDeltaSigma: Text generation
#-------------------------------------------------------------------------------

const ORDINAL_STRING = [
	"Zeroth", "First", "Second", "Third", "Fourth", "Fifth",
	"Sixth", "Seventh", "Eighth", "Nineth", "Tenth",
	"Twelfth", "Thirteenth",
]

function ds_orderString(n::Int, capitalize::Bool=false)
	idx = n+1
	if idx in keys(ORDINAL_STRING)
		s = ORDINAL_STRING[idx]
		capitalize || (s = lowercase(s))
	else
		s = "$(n)th"
	end
	return s
end


function str_modulatortype(dsm::AbstractDSM)
	typestr = (0==dsm.f0) ? "Lowpass" : "Bandpass"
	if isquadrature(dsm)
		typestr = "Quadrature " * typestr
	end
	return ds_orderString(dsm.order, true) * "-Order " * typestr
end

#Last line
