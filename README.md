<!-- Reference-style links to make tables & lists more readable -->
[Gallery]: <https://github.com/ma-laforge/FileRepo/tree/master/SignalProcessing/sampleplots/README.md>
[CMDimDataJL]: <https://github.com/ma-laforge/CMDimData.jl>
[InspectDRJL]: <https://github.com/ma-laforge/InspectDR.jl>


# RSDeltaSigmaPort.jl: Port of Richard Schreier's Delta Sigma Toolbox&sup1;

&sup1;Richard Schreier (2021). Delta Sigma Toolbox (<https://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox>), MATLAB Central File Exchange. Retrieved March 20, 2021.

[![Build Status](https://travis-ci.org/ma-laforge/RSDeltaSigmaPort.jl.svg?branch=master)](https://travis-ci.org/ma-laforge/RSDeltaSigmaPort.jl)

***INITIAL STAGES OF PORT***: Sorry. Not yet ready for use.

## Table of contents

 1. [Description](#Description)
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

Extensive compatibility testing of `RSDeltaSigmaPort.jl` has not been performed.  The module has been tested using the following environment(s):

- Linux / Julia-1.5.3

## Disclaimer

 - ***INITIAL STAGES OF PORT***: Sorry. Not yet ready for use.
 - The `RSDeltaSigmaPort.jl` module is not yet mature.  Expect significant changes.
