# RSDeltaSigma: TODO

## Migrate to use ControlSystems.jl directly?
Deprecate legacy functions/names to use slightly altered API from ControlSystems.jl?
```julia
#See:
export _ss, _zpk, _zpkdata
export _zp2ss, _zp2tf, _freqz, _tf, _minreal, _impulse
#etc.
```

## Deprecate `array_round()`
 - Figure out if orignial code rounded indicies to nearest integer, or implicitly applied `floor()`, etc.
 - Start using `div(num, den)` more if possible (should be faster than `round`/`floor`/`ceil`).

## Add regression tests.
Very little code is tested at the moment.

## List other known elements needing attention HERE.

