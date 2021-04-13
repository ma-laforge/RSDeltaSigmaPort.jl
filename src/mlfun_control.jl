#RSDeltaSigmaPort: Implement missing core functionality using ControlSystems.jl
#-------------------------------------------------------------------------------

#=TODO (Goal)
Want to not depend on ControlSystems.jl (slows down compile time significantly).
=#

function _zp2ss(z, p, k)
	H = ControlSystems.zpk(z, p, k)
	return ControlSystems.ss(H)
end

function _zp2tf(z, p, k)
	H = ControlSystems.zpk(z, p, k)
	_num = ControlSystems.num(H)
	_den = ControlSystems.den(H)
	return (_num, _den)
end

function _freqz(num, den, w)
	H = ControlSystems.tf(num, den, 1)
	return ControlSystems.freqresp(H, w)
end

_tf(num, den, ts) = ControlSystems.tf(num, den, ts)

function _impulse(sys, Tfinal::Real)
	y, t, x = ControlSystems.impulse(sys, Float64(Tfinal))
	return y
end

#Last line
