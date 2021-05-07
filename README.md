<!-- Reference-style links to make tables & lists more readable -->
[Gallery]: <https://github.com/ma-laforge/FileRepo/blob/master/RSDeltaSigmaPort/notebook>
[CMDimDataJL]: <https://github.com/ma-laforge/CMDimData.jl>
[InspectDRJL]: <https://github.com/ma-laforge/InspectDR.jl>


# RSDeltaSigmaPort.jl: Port of Richard Schreier's Delta Sigma Toolbox&sup1;
**Galleries:** [:art: Sample notebooks (w/outputs)][Gallery]

&sup1;Richard Schreier (2021). Delta Sigma Toolbox (<https://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox>), MATLAB Central File Exchange. Retrieved March 20, 2021.

[![Build Status](https://github.com/ma-laforge/RSDeltaSigmaPort.jl/workflows/CI/badge.svg)](https://github.com/ma-laforge/RSDeltaSigmaPort.jl/actions?query=workflow%3ACI)

### :warning: Progress report
***INTERMEDIATE STAGE OF PORT***: A significant portion of the Delta Sigma toolbox has been ported.

The following high-level functionnality has (at least partially) been ported:
 - `simulateDSM`, `simulateMS`, `simulateSNR`, `simulateHBF`
 - `synthesizeNTF`, `realizeNTF`, `realizeNTF_ct`
 - `calculateSNR`, `peakSNR`, `predictSNR`
 - `calculateTF`, `evalTF`, `evalTFP`
 - `stuffABCD`, `scaleABCD`, `mapABCD`, `partitionABCD`
 - `mapCtoD`, `mapQtoR`
 - `exampleHBF`
 - `pulse`, `impL1`
 - `lollipop`, `logsmooth`
 - `documentNTF`, `plotExampleSpectrum`

And demos:
 - `dsdemo1`, ..., `dsdemo6`, `dsexample1`, `dsexample2`, `demoLPandBP`

## Table of contents

 1. [Description](#Description)
 1. Sample usage
    1. [IJupyter notebooks (`notebook/`)](notebook/)
    1. [IJupyter notebooks (with output) &#x21AA;][Gallery]
    1. [Sample directory w/plain `.jl` files (`sample/`)](sample/)
 1. [Plotting](#Plotting)
 1. [Installation](#Installation)
 1. [Running sample scripts](#SampleScripts)
 1. [Julia tips](doc/juliatips.md)
    1. [Useful functions](doc/juliatips.md#FunctionLibraries)
    1. [`linspace()` & `logspace()`](doc/juliatips.md#LinLogSpace)
 1. [Known limitations](#KnownLimitations)
    1. [TODO](doc/todo.md)


<a name="Description"></a>
## Description
As its name suggests, `RSDeltaSigmaPort.jl` is a Julia port of Richard Schreier's Delta Sigma Toolbox.

### Module name
Note that this module is not named something like `DeltaSigmaModulators.jl`, thus allowing someone else to appropriate the package name later on. Hopefully, `RSDeltaSigmaPort` will eventually be superseded by a more generically named package that is better integrated with Julia's ecosystem than this simple port.

### Design decisions
This module tries to remain true to the original Delta Sigma Toolbox while conforming to some Julia patterns, including:
 - Multiple dispatch (make function calls simpler to write).
 - ***Not*** writing each function definition in its own, separate file.
 - Using keyword arguments when deemed appropriate.
 - Returning `NamedTuple`s instead of simple arrays when multiple values are returned.
 - ...

Progressively replacing modulator parameters in function calls with `RealDSM` and `QuadratureDSM` objects:
 - Simplifies function interface for user.
 - Centralizes defaults for parameter values on construction of `RealDSM` and `QuadratureDSM`.
 - Looking to keep "original" function interface (with individual modulator parameters) available for accustomed users.
 - Looking to remove default values from said interface to avoid unexpected bugs from inconsistent defaults.
 - Might change with time (not sure if certain parameters, like `opt`, should migrate to a NTF structure or something).

<a name="Plotting"></a>
## Plotting
`RSDeltaSigmaPort.jl` uses [CMDimData.jl/EasyPlot][CMDimDataJL] to handle plotting.
For examples on how to generate ***new/customized*** plots, see the built-in
functions found in `plot_*.jl` files in the source directory:
 - [`src/`](src/)

<a name="Installation"></a>
## Installation
The `RSDeltaSigmaPort.jl` toolbox is written using the Julia programming
language. Unless you already have Julia installed, you will need to first
install the base language. Simply download \& install the most recent version
of Julia from Julia's official "downloads" page.

**Julia's official "downloads" page:**
 - <https://julialang.org/downloads/>

Step 2 is to install the `RSDeltaSigmaPort.jl` package itself. Since
`RSDeltaSigmaPort.jl` is registered with Julia's **General** registry, you can
automatically download & install it from Julia's built-in package manager.
Simply launch Julia, and run the following from the command prompt:

```julia-repl
julia> ]
pkg> add RSDeltaSigmaPort
```

<a name="SampleScripts"></a>
## Running sample scripts
Sample scripts in the `sample/` subdirectory can be run using `include()`.

For convenience, the `@runsample` macro automatically locates the script path
and executes `include()` for you:

```julia-repl
julia> using RSDeltaSigmaPort #Will take a while to load, compile, etc...
julia> import RSDeltaSigmaPort: @runsample

julia> @runsample("dsdemo1.jl")
julia> @runsample("dsdemo2.jl")
julia> @runsample("dsdemo3.jl")
julia> @runsample("dsdemo4_audio.jl")
julia> @runsample("dsdemo5.jl")
julia> @runsample("dsdemo6.jl")
julia> @runsample("dsexample1.jl")
julia> @runsample("dsexample2.jl")
julia> @runsample("demoLPandBP.jl")
```

<a name="KnownLimitations"></a>
## Known limitations
Functions that are not supported:
 - `printmif()`

### [TODO](doc/todo.md)

### Compatibility

Extensive compatibility testing of `RSDeltaSigmaPort.jl` has not been performed.
The module has been tested using the following environment(s):

- Linux / Julia-1.6.0

## :warning: Disclaimer

 - ***INTERMEDIATE STAGE OF PORT***: A significant portion of the Delta Sigma toolbox has been ported.
 - Jupyter [notebooks](notebook/) might be slightly broken/out of date. If so,
   see their counterparts in the [`sample/`](sample/) directory for a more
   regularly maintained example.
 - The `RSDeltaSigmaPort.jl` module is not yet mature.  Expect significant changes.
