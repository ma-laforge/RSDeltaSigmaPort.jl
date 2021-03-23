# RSDeltaSigmaPort.jl: Port of Richard Schreier's Delta Sigma Toolbox&sup1;

&sup1;Richard Schreier (2021). Delta Sigma Toolbox (<https://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox>), MATLAB Central File Exchange. Retrieved March 20, 2021.

***INITIAL STAGES OF PORT***: Sorry. Not yet ready for use.

This module tries to remain true to the original Delta Sigma Toolbox while conforming to some Julia patterns, including:
 - ***Not*** writing each function definition in its own, separate file.
 - Using keyword arguments when deemed appropriate.
 - ...

### Module name
Note that this module is not named `DeltaSigma.jl`, thus allowing someone else to appropriate the name later on. The intent is to eventually have an "official" `DeltaSigma.jl` module that is better integrated with Julia's ecosystem than this simple port.

<a name="KnownLimitations"></a>
## Known limitations

### Compatibility

Extensive compatibility testing of `RSDeltaSigmaPort.jl` has not been performed.  The module has been tested using the following environment(s):

- Linux / Julia-1.5.3

## Disclaimer

 - ***INITIAL STAGES OF PORT***: Sorry. Not yet ready for use.
 - The `RSDeltaSigmaPort.jl` module is not yet mature.  Expect significant changes.
