# Hohmann Transfer Visualization
This is a script written in the [Julia](https://julialang.org/) programming language for visualizing Hohmann transfers.

![](/hohmann_transfer_example.mp4)

This script uses the [Makie](https://github.com/JuliaPlots/Makie.jl) package to handle the visualization.

To run the visualization make sure you have the Julia programming language installed along with the following packages:
- CSV
- GeometryBasics
- LinearAlgebra
- GLMakie
- FileIO
- ColorSchemes
- DataFrames
- OrdinaryDiffEq
- JLD2

To run the simulation, first clone the repository and change into the `src` directory:
```bash
git clone https://github.com/crispy-landslide/Hohmann-Transfer-Visualization.git
cd Hohmann-Transfer-Visualization/src
```

Then, start a Julia session
```bash
julia
```
Inside your Julia session run the following:
```julia
include("main.jl")
```

It will take some time (over 60 seconds) to run the first time. After the simulation is run once, it will be much faster (less than 1 second) inside the same session.




