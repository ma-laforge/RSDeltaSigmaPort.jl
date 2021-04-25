# RSDeltaSigma: Julia tips

 1. [Useful functions](#FunctionLibraries)
 1. [`linspace()` & `logspace()`](#LinLogSpace)

<a name="FunctionLibraries"></a>
## Useful functions
Not all functions you might want to use are directly available in julia.
In many cases, you need to make functions available by importing them from
different libraries:

```julia
import Random
import Statistics: mean
import LinearAlgebra
import LinearAlgebra: norm, diagm, eigen
import SpecialFunctions: erfinv
import FFTW: fft, fftshift
import DSP: conv
import Polynomials: roots, Polynomial
import Printf: @sprintf
```

`RSDeltaSigma` does not automatically export these symbols/functions (make
the functions available to the user). This choice was deliberatly made to
avoid symbol collisions in cases where the user needs to manually import
said symbols to get other portions of their code working.

<a name="LinLogSpace"></a>
## `linspace()` & `logspace()`
Functions `linspace()` & `logspace()` don't exist in Julia 1.0+. Instead, you should use `range()`:

 1. `range(start, stop=stop, length=n)`: Constructs a `StepRangeLen` object.
     - `collect(range(0, stop=100, length=101))`: Generates a `Vector{Float64}` (not `Vector{Int}`).
 1. `10 .^ range(start, stop=stop, length=n)`: Generates a `Vector{Float64}` object.
     - `10 .^ range(-10, stop=10, length=101)`: Log-spaced values &isin; [1.0e-10, 1.0e10].


