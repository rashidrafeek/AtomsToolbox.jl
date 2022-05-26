module AtomsToolbox

using Distances: Distances, pairwise, euclidean, peuclidean,
                 Euclidean, PeriodicEuclidean
using AtomsBase: AbstractSystem,
                 FastSystem,
                 FlexibleSystem,
                 position,
                 atomic_number,
                 atomic_symbol,
                 atomic_mass,
                 bounding_box,
                 boundary_conditions,
                 periodicity
using Unitful: Unitful, ustrip, @u_str, unit
using Graphs: SimpleGraph, connected_components
using LinearAlgebra: LinearAlgebra, ⋅, det, norm
using StaticArrays: StaticMatrix, Size

# Constants
export covalent_radii
include("constants.jl")

# Integration with AtomsBase
include("atomsbase.jl")

# Functions to get data from `AbstractSystem`s
export box_lengths,
       getdistance,
       getdistancematrix,
       getangle,
       getvolume,
       cell_angles,
       cell_lengths_and_angles,
       natural_cutoffs,
       getconnectivitymatrix,
       getconnectedcomponents,
       transformpositions,
       wrap
include("getters.jl")

# Functions to transform `AbstractSystem`s
export transformpositions, wrap
include("transform.jl")

# Fuction to build new `AbstractSystem`s
include("build.jl")

end
