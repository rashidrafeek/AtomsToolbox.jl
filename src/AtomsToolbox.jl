module AtomsToolbox

using Distances: Distances, pairwise, euclidean, peuclidean,
                 Euclidean, PeriodicEuclidean
using AtomsBase: AtomsBase,
                 AbstractSystem,
                 FastSystem,
                 FlexibleSystem,
                 Atom,
                 AtomView,
                 position,
                 atomic_number,
                 atomic_symbol,
                 atomic_mass,
                 cell_vectors,
                 periodicity,
                 species
using Unitful: Unitful, ustrip, @u_str, unit
import Graphs
using LinearAlgebra: LinearAlgebra, ⋅, det, norm, lu, ×
using StaticArrays: StaticMatrix, Size
import Base: angle, sort # To extend for AbstractSystem 

# Some type piracy to deal with other packages
# Should be upstreamed or removed eventually
include("piracy.jl")

# Constants
include("constants.jl")

# Integration with AtomsBase
include("atomsbase.jl")

# Functions to get data from `AbstractSystem`s
export distance,
       distance_matrix,
       angle,
       dihedral,
       volume,
       cell_lengths,
       cell_angles,
       cell_parameters,
       scaled_position
include("getters.jl")

# Functions to transform `AbstractSystem`s
export wrap,
       supercell
include("transform.jl")

# Function to build new `AbstractSystem`s
include("build.jl")

end
