#RSDeltaSigmaPort: Implement missing core functionality using ControlSystems.jl
#-------------------------------------------------------------------------------

#=TODO (Goal)
Want to not depend on ControlSystems.jl (slows down compile time significantly).

UPDATE: Much better in Julia 1.6 - especially using Revise.jl.
Maybe ok to depend on ControlSystems.jl.
=#

_c2d(sys, Ts, method=:zoh) = ControlSystems.c2d(sys, Ts, method)

_ss(A,B,C,D) = ControlSystems.ss(A,B,C,D) #Ts=0.0 default: continuous-time system
_ss(A,B,C,D,Ts) = ControlSystems.ss(A,B,C,D,Ts)

function _zpk(sys::ControlSystems.StateSpace)
	H=ControlSystems.zpk(sys)
	z, p, k = ControlSystems.zpkdata(H)

	#Cannot access sys.Ts for continuous time:
	#(use sys because Ts does not get ported over in zpk() for some reason)
	Ts = isdiscrete(sys) ? sys.Ts : 0.0
	return _zpk(z, p, k, Ts)
end

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

_minreal(sys, tol) = ControlSystems.minreal(sys, tol)

function _impulse(sys, t)
	y, t, x = ControlSystems.impulse(sys, t)
	return (y, t, x)
end

function _lsim(sys, u)
	if !isdiscrete(sys)
		throw("_lsim: must specify t-vector")
	end
	Ts = sys.Ts
	t = range(0, step=Ts, length=length(u))
	y, t, x = ControlSystems.lsim(sys, u, t)
	return (y, t, x)
end


#Last line
