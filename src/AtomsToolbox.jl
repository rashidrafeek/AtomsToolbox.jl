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
import Graphs
using LinearAlgebra: LinearAlgebra, ⋅, det, norm
using StaticArrays: StaticMatrix, Size
import Base: angle # To extend for AbstractSystem 

# Constants
include("constants.jl")

# Integration with AtomsBase
include("atomsbase.jl")

# Functions to get data from `AbstractSystem`s
export distance,
       distance_matrix,
       angle,
       volume,
       cell_lengths,
       cell_angles,
       cell_parameters,
       scaled_position
include("getters.jl")

# Functions to transform `AbstractSystem`s
export wrap
include("transform.jl")

# Fuction to build new `AbstractSystem`s
include("build.jl")

end
