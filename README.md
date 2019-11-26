![Ju-SFEM](/docs/images/logo.png)

<div align="center">
<h1><b><i>juSFEM</i></b> User Manual</h1>

[![](https://img.shields.io/badge/Lang-Julia_v1.1.1-brightgreen.svg)](https://julialang.org/)&nbsp;
[![](https://badgen.net/badge/icon/Juno?icon=atom&label)](https://junolab.org/)&nbsp;
[![](https://img.shields.io/badge/Platform-windows_10_x64-brightgreen.svg)](https://www.microsoft.com/en-us/windows)&nbsp;
[![](https://img.shields.io/badge/Licence-MIT-blue.svg)](https://github.com/Elliothuo/Ju-SFEM/blob/master/LICENSE)

Numerical simulation package written in julia

</div>

> This is the development version of the source code, the function is not complete, only to achieve the preliminary work.
`2019.06.14`

* 1.Introduction
* 2.Structure of the Package juSFEM
* 3.The Developed Package juSFEM
* 4.Validation and Evaluation of juSFEM
* 5.Computational Correctness and Speedup
* 6.How to Use

# 1. Introduction

To balance the program execution **efficiency** and the **ease of
implementing** the program, we design and implement an open-source
package of parallel S-FEM for elastic problems by using the Julia
language. We term our package as juSFEM.

# 2. Structure of the Package juSFEM

<div align="center"> <img
src="/docs/images/figure1.png"
alt="Figure 1" width="500" height="300"></img>
</div>

> Structure of the package juSFEM

# 3. The Developed Package juSFEM

- Pre-processing

<div align="center"> <img
src="/docs/images/figure2.png"
alt="Figure 2" width="500" height="300"></img>
</div>

>  Illustraction of processing .face file (Take two adjacent tetrahedrons as an example.)

- Solver

<div align="center"> <img
src="/docs/images/figure3.png"
alt="Figure 3" width="500" height="300"></img>
</div>

> Process of the Solver

- Post-processing

<div align="center"> <img
src="/docs/images/figure4.png"
alt="Figure 4" width="350" height="200"></img>
</div>

> Plotting with the ParaView

<div align="center"> <img
src="/docs/images/figure5.gif"
alt="Figure 5" width="680" height="450"></img>
</div>

> Plotting of the displacement with the PlotlyJS

# 4. Validation and Evaluation of juSFEM

<div align="center">

| Specifications | Details |
| :----: | :----: |
| OS | Windows 10 Professional |
| CPU | Intel Xeon Gold 5118 CPU |
| CPU Frequency | 2.30 (GHz) |
| CPU Cores  | 24 |
| Julia Version  | Julia 1.1.1 |

</div>

## 4.1 Verification of the Accuracy of juSFEM

In this example, the size of the 342
cantilever beam is 1×0.2×0.2(m). The upper face of the cantilever
beam is subjected to a uniform pressure p = 12500(N/m^2) in the
direction of Z.The elastic modulus E = 2E8(N/m^2) and Poisson’s
ratio v = 0.3. The model consists of 458 nodes and 1,576
tetrahedral elements.

<div align="center"> <img
src="/docs/images/figure6.png"
alt="Figure 6" width="450" height="200"></img>
</div>

<div align="center"> <img
src="/docs/images/figure7.png"
alt="Figure 7" width="450" height="200"></img>
</div>

>  (a) A 3D cantilever beam subjected to a uniform pressure on the top face and (b) a mesh
with four-node tetrahedral elements.

<div align="center"> <img
src="/docs/images/figure8.png"
alt="Figure 8" width="450" height="440"></img>
</div>

<div align="center"> <img
src="/docs/images/figure9.png"
alt="Figure 9" width="650" height="300"></img>
</div>

>  Comparison of the Cantilever beam deflection curves

## 4.2 Evaluation of the Efficiency of juSFEM

<div align="center">

| Mesh model (T4) | Number of Nodes | Number of Elements |
| :----: | :----: | :----: |
| 1 | 94,066 | 503,963 |
2 | 191,132 | 1,037,691 |
3 | 275,655 | 1,498,126 |
4 | 382,375 | 2,076,930 |

</div>

<div align="center"> <img
src="/docs/images/figure10.png"
alt="Figure 10" width="450" height="320"></img>
</div>

> Comparison of the computational time for assembling the global stiffness matrix

<div align="center"> <img
src="/docs/images/figure11.png"
alt="Figure 11" width="450" height="320"></img>
</div>

> Comparison of the overall computational time of the parallel juSFEM on multi-core CPU,
the serial juSFEM on single-core CPU, and the software ABAQUS

# 5. Computational Correctness and Speedup

<div align="center">

| Method | 0.2m | 0.4m | 0.6m | 0.8m | 1.0m |
| :----: | :----: | :----: | :----: | :----: | :----: |
| FEM-T4 | -8.20E-04 | -2.61E-03 | -4.99E-03 | -7.62E-03 | -1.03E-02 |
| juSFEM-T4 | -8.66E-04 | -2.82E-03 | -5.41E-03 | -8.28E-03 | -1.12E-02 |
| FEM-T10 | -9.54E-04 | -3.08E-03 | -5.86E-03 | -8.92E-03 | -1.20E-02 |

</div>

<div align="center"> <img
src="/docs/images/figure12.png"
alt="Figure 12" width="450" height="320"></img>
</div>

> The speedup of assembling global stiffness matrix in parallel

<div align="center"> <img
src="/docs/images/figure13.png"
alt="Figure 13" width="600" height="200"></img>
</div>

> The proportion of the time spent in assembling stiffness matrix and solving the equation
(a. single-core b. multi-core)

# 6. How to Use

## 6.1 multi-core juSFEM

`src`->`multi-core juSFEM`->`SFEM.jl`
```julia
using Distributed
addprocs(your computer CPU processors number)
...
include("path/to/assembleout.jl")
include("path/to/assemblein.jl")
```

`src`->`multi-core juSFEM`->`assemblein.jl`
```julia
include("path/to/area.jl")
include("path/to/vectorin.jl")
include("path/to/volume.jl")
```

`src`->`multi-core juSFEM`->`assembleout.jl`
```julia
include("path/to/area.jl")
include("path/to/vectorout.jl")
include("path/to/volume.jl")
```

`src`->`multi-core juSFEM`->`main.jl`
```julia
A = "path/to/4.cen"
B = "path/to/4.ele"
C = "path/to/4.facein"
D = "path/to/4.faceout"
E = "path/to/4.node"
```

## 6.2 single-core juSFEM

`src`->`single-core juSFEM`->`SFEM.jl`
```julia
include("path/to/area.jl")
include("path/to/vectorout.jl")
include("path/to/vectorin.jl")
include("path/to/volume.jl")
include("path/to/assemblehex.jl")
include("path/to/assembletet.jl")
```

`src`->`single-core juSFEM`->`main.jl`
```julia
A = "path/to/test.cen"
B = "path/to/test.ele"
C = "path/to/test.facein"
D = "path/to/test.faceout"
E = "path/to/test.node"
F = "path/to/test.face"
```

> The test files can be downloaded [here](http://www.google.com/).

# To-do List

- [x] S-FEM 2D (ES-FEM)
- [x] S-FEM 3D (FS-FEM)
- [x] The elastic problems on multi-core processors
- [ ] The plastic problems on multi-core processors
- [ ] GPU parallel in juSFEM
- [ ] Practical application
- [ ] Post-processing methods
- [ ] Register Julia package
