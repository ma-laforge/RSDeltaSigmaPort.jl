#RSDeltaSigmaPort: Implement missing core functionality using ControlSystems.jl
#-------------------------------------------------------------------------------

#=TODO (Goal)
Want to not depend on ControlSystems.jl (slows down compile time significantly).

UPDATE: Much better in Julia 1.6 - especially using Revise.jl.
Maybe ok to depend on ControlSystems.jl.
=#

_ss(A,B,C,D,Ts=1) = ControlSystems.ss(A,B,C,D,Ts)

function _zpk(sys::ControlSystems.StateSpace)
	H=ControlSystems.zpk(sys)
	z, p, k = ControlSystems.zpkdata(H)
	Ts = sys.Ts #Does not get ported over in zpk() for some reason
	#VERIFYME Not sure why data is on 2nd index
	#when StateSpace constructed with _ss() function
	z=z[2][:]; p=p[2][:]; k=k[2]
#	map(display, [z, p, k])
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

function _impulse(sys, Tfinal::Real)
	y, t, x = ControlSystems.impulse(sys, Float64(Tfinal))
	return y
end

function _minreal(tf::ZPKData, tol)
	tf = ControlSystems.zpk(tf.z, tf.p, tf.k, tf.Ts)
	mrtf = ControlSystems.minreal(tf, tol)
	z, p, k = ControlSystems.zpkdata(mrtf)
	z=z[1][:]; p=p[1][:]; k=k[1]
	Ts = mrtf.Ts
#	map(display, [z, p, k, Ts])
	return _zpk(z, p, k, Ts)
end



#Last line
