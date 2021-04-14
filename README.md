<!-- Reference-style links to make tables & lists more readable -->
[Gallery]: <https://github.com/ma-laforge/FileRepo/blob/master/RSDeltaSigmaPort/notebook>
[CMDimDataJL]: <https://github.com/ma-laforge/CMDimData.jl>
[InspectDRJL]: <https://github.com/ma-laforge/InspectDR.jl>


# RSDeltaSigmaPort.jl: Port of Richard Schreier's Delta Sigma Toolbox&sup1;
**Galleries:** [:art: Sample notebooks (w/outputs)][Gallery]

&sup1;Richard Schreier (2021). Delta Sigma Toolbox (<https://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox>), MATLAB Central File Exchange. Retrieved March 20, 2021.

[![Build Status](https://github.com/ma-laforge/RSDeltaSigmaPort.jl/workflows/CI/badge.svg)](https://github.com/ma-laforge/RSDeltaSigmaPort.jl/actions?query=workflow%3ACI)

### :warning: Progress report
***INITIAL STAGES OF PORT***: A limited portion of the Delta Sigma toolbox has been ported.

The following high-level functionnality has (at least partially) been ported:
 - `evalTF`
 - `calculateSNR`, `peakSNR`, `predictSNR`
 - `synthesizeNTF`
 - `simulateSNR`, `simulateDSM`

## Table of contents

 1. [Description](#Description)
 1. Sample usage
    1. [IJupyter notebooks (`notebook/`)](notebook/)
    1. [IJupyter notebooks (with output) &#x21AA;][Gallery]
    1. [Sample directory w/plain `.jl` files (`sample/`)](sample/)
 1. [Plotting](#Plotting)
 1. [Installation](#Installation)
 1. [Known limitations](#KnownLimitations)


<a name="Description"></a>
## Description
As its name suggests, `RSDeltaSigmaPort.jl` is a Julia port of Richard Schreier's Delta Sigma Toolbox.

### Module name
Note that this module is not named `DeltaSigma.jl`, thus allowing someone else to appropriate the name later on. The intent is to eventually have an "official" `DeltaSigma.jl` module that is better integrated with Julia's ecosystem than this simple port.

### Design decisions
This module tries to remain true to the original Delta Sigma Toolbox while conforming to some Julia patterns, including:
 - ***Not*** writing each function definition in its own, separate file.
 - Using keyword arguments when deemed appropriate.
 - ...

<a name="Plotting"></a>
## Plotting
`RSDeltaSigmaPort.jl` uses [CMDimData.jl/EasyPlot][CMDimDataJL] to handle plotting.
For examples on how to generate ***new/customized*** plots, see the built-in
functions found in:
 - [plotgen.jl](src/plotgen.jl).

<a name="Installation"></a>
## Installation
`RSDeltaSigmaPort.jl` is not yet registered with Julia's **General** registry.
It can nonetheless be installed using the library's URL from Julia's built-in package manager:

```julia-repl
julia> ]
pkg> add https://github.com/ma-laforge/RSDeltaSigmaPort.jl
```

<a name="KnownLimitations"></a>
## Known limitations

### Compatibility

Extensive compatibility testing of `RSDeltaSigmaPort.jl` has not been performed.
The module has been tested using the following environment(s):

- Linux / Julia-1.6.0

## :warning: Disclaimer

 - ***INITIAL STAGES OF PORT***: A limited portion of the Delta Sigma toolbox has been ported.
 - Jupyter [notebooks](notebook/) might be slightly broken/out of date. If so,
   see their counterparts in the [`sample/`](sample/) directory for a more
   regularly maintained example.
 - The `RSDeltaSigmaPort.jl` module is not yet mature.  Expect significant changes.
