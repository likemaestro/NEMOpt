# NEMOpt: Nonlinear Elastic Microstructural Topology Optimization

NEMOpt is a Julia-based research codebase for topology optimization of periodic nonlinear elastic microstructures, with a focus on designing metamaterials that exhibit targeted macroscopic behavior (for example, negative Poisson's ratio / auxetic response).

The repository provides two workflows:

- Full-domain optimization (`Full_Domain`)
- Domain-reduction optimization (`Domain_Reduction`)

Both workflows use nonlinear finite element analysis, density filtering/projection, homogenization, and MMA-based constrained optimization.

## Highlights

- Nonlinear microstructural analysis under prescribed macroscopic deformation gradient.
- Helmholtz filtering with Heaviside projection for robust topology updates.
- Homogenized property tracking during optimization:
	- Poisson's ratio (`v_12`)
	- Bulk modulus (`K`)
	- Shear modulus (`G`)
- MMA optimization (`MMA02` / `MMA87`) via `Nonconvex.jl`.
- Iteration-wise topology export and convergence plots.
- Animated topology evolution export (`Optimum-Topology.mp4`).

## Repository Structure

```text
NEMOpt/
|-- Full_Domain/
|   |-- NEMOpt_main.jl         # Entry point: full-domain run
|   |-- NEMOpt.jl              # Full-domain implementation
|   `-- Initial_Designs/       # Input initial topology images
|-- Domain_Reduction/
|   |-- NEMOpt-DR_main.jl      # Entry point: domain-reduction run
|   |-- NEMOpt-DR.jl           # Domain-reduction implementation
|   `-- Initial_Designs/       # Input initial topology images
`-- README.md
```

## Requirements

- Julia 1.9+ (recommended)
- An `ffmpeg` executable available to Julia/FFMPEG.jl for video export
- Julia packages used by this project:
	- `Nonconvex`, `ChainRulesCore`, `Dates`, `Random`, `Distributions`
	- `SparseArrays`, `LinearAlgebra`, `Plots`, `Colors`, `Images`
	- `JLD2`, `DataFrames`, `Printf`, `DelimitedFiles`, `Tullio`
	- `FFMPEG`, `PlotlyJS`

## Setup

From the repository root, install dependencies in Julia:

```julia
using Pkg
Pkg.add([
		"Nonconvex", "ChainRulesCore", "Distributions", "Plots", "Colors", "Images",
		"JLD2", "DataFrames", "Tullio", "FFMPEG", "PlotlyJS"
])
Pkg.add("NonconvexMMA")
```

Note: Standard-library modules such as `Dates`, `Random`, `SparseArrays`, `LinearAlgebra`, `Printf`, and `DelimitedFiles` do not need installation.

## Quick Start

### 1) Full-domain workflow

```bash
cd Full_Domain
julia NEMOpt_main.jl
```

### 2) Domain-reduction workflow

```bash
cd Domain_Reduction
julia NEMOpt-DR_main.jl
```

Each run creates a timestamped output folder in the active working directory.

## Input Configuration

Main tunable parameters are defined near the top of each entry script:

- `Full_Domain/NEMOpt_main.jl`
- `Domain_Reduction/NEMOpt-DR_main.jl`

Typical parameters:

- Geometry and discretization: `Lx`, `Ly`, `Nelx`, `Nely`
- Material: `E`, `v`
- Loading: `F_M`
- Filtering/projection: `r`, `p`, `epsilon`, `beta`, `omega`, `beta2`, `omega2`
- Optimization controls: `N_ITER`, `V_MAX`
- Initialization choice (exactly one):
	- `INITIAL_DESIGN_NAME = "...png"`
	- `RANDOM_SEED = ...`

For image-based initialization, place topology images in the corresponding `Initial_Designs` folder.

## Outputs

For each run, NEMOpt generates a timestamped workspace including:

- `Parameters.txt`, `Description.txt`
- `Topologies/` (iteration snapshots)
- `Graphs/` (objective and effective-property history plots)
- `Homogenized_Stiffnesses/` and `Homogenized_Stresses/` per iteration
- `Elapsed_Time.txt`
- `Optimum-Topology.mp4`

## Citation

If you use NEMOpt in academic work, please cite: Guven, M. and Ozdemir, I., "Topology Optimization of Nonlinear Elastic Micro-structures by a Domain Reduction Technique"