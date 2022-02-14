module AtomsToolbox

using Distances: pairwise,
                 euclidean,
                 peuclidean,
                 Euclidean,
                 PeriodicEuclidean
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
using Unitful: ustrip,
               @u_str
using Graphs: SimpleGraph, 
              connected_components

export covalent_radii, 
       box_lengths,
       getdistance,
       getdistancematrix,
       natural_cutoffs, 
       getconnectivitymatrix,
       getconnectedcomponents,
       transformpositions,
       wrap

const covalent_radii = [
                        0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57,
                        0.58, 1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06,
                        2.03, 1.76, 1.7,  1.6,  1.53, 1.39, 1.39, 1.32, 1.26,
                        1.24, 1.32, 1.22, 1.22, 1.2,  1.19, 1.2,  1.2,  1.16,
                        2.2,  1.95, 1.9,  1.75, 1.64, 1.54, 1.47, 1.46, 1.42,
                        1.39, 1.45, 1.44, 1.42, 1.39, 1.39, 1.38, 1.39, 1.4,
                        2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98,
                        1.96, 1.94, 1.92, 1.92, 1.89, 1.9,  1.87, 1.87, 1.75,
                        1.7,  1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, 1.45,
                        1.46, 1.48, 1.4,  1.5,  1.5,  2.6,  2.21, 2.15, 2.06,
                        2.0,  1.96, 1.9,  1.87, 1.8,  1.69, 0.2,  0.2,  0.2,
                        0.2,  0.2,  0.2,  0.2,  0.2,  0.2,  0.2,  0.2,  0.2,
                        0.2,  0.2,  0.2,  0.2,  0.2,  0.2,  0.2,  0.2,  0.2,
                        0.2
                       ]

# function Base.getindex(system::AbstractSystem, range::AbstractVector)
#     positions = position(system)[range]
#     cell = bounding_box(system)
#     at = species(system)[range]
# 
#     typeof(system)(cell, positions, at)
# end

"""
    _getconstructor(s::AbstractSystem)

Get the constructor to construct the concrete type of `s` which is a subtype of
AbstractSystem
"""
_getconstructor(s::AbstractSystem) = error("""
    No constructor defined. Please define _getconstructor(::$(Base.typename(typeof(s)).name)).
    """)
_getconstructor(::FastSystem) = FastSystem
_getconstructor(::FlexibleSystem) = FlexibleSystem

function Base.getindex(system::AbstractSystem, range::AbstractVector)
    particles = collect(system)
    _getconstructor(system)(
                            particles[range], bounding_box(system), boundary_conditions(system)
                           )
    # _getconstructor(system)(
    #     bounding_box(system), boundary_conditions(system), 
    #     position(system)[range], atomic_symbol(system)[range],
    #     atomic_number(system)[range], atomic_mass(system)[range])
end

box_lengths(sys::AbstractSystem) = getindex.(bounding_box(sys),[1,2,3])

"""
    getdistance(system::AbstractSystem, at1, at2)

Get the distance between atoms with indices `at1` and `at2`.
"""
function getdistance(system::AbstractSystem, at1, at2)
    pos1 = position(system, at1)
    pos2 = position(system, at2)
    if all(periodicity(system))
        cell = box_lengths(system)
        dist = peuclidean(pos1, pos2, cell)
    else
        dist = euclidean(pos1, pos2)
    end

    return dist
end

"""
    getdistancematrix(system::AbstractSystem)

Get the distance matrix between all atoms as an NxN matrix where N is the
number of atoms in the given `system`.
"""
function getdistancematrix(system::AbstractSystem)
    if all(periodicity(system))
        pairwise(
            PeriodicEuclidean(ustrip.(box_lengths(system))),
            ustrip.(position(system))
        )
    else
        pairwise(Euclidean(), ustrip.(position(system)))
    end
end

"""
    natural_cutoffs(system::AbstractSystem)

Get the covalent radii of the atoms in the given `system`. Currently uses values
obtained from ASE.
"""
function natural_cutoffs(system::AbstractSystem)
    radii = covalent_radii[atomic_number.(system)]
    return radii
end

"""
    getconnectivitymatrix(
        system;                                                                 
        nlcutoff = natural_cutoffs(system).+0.2,
        retdistmat = false
    )

Get the connectivity matrix of the given `system`. The cutoff is specified
using `nlcutoff`. To return the distance matrix, set `retdistmat` as `true`.
"""
function getconnectivitymatrix(
        system;
        nlcutoff = natural_cutoffs(system).+0.2,
        retdistmat = false
    )   
    distmat = getdistancematrix(system)
    
    connmat = (distmat .- nlcutoff') .< nlcutoff
    
    if retdistmat
        return connmat, distmat
    else
        return connmat
    end
end

"""
    getconnectedcomponents(
        system, dists=true;
        nlcutoff = natural_cutoffs(str).+0.2,
        retconnmat = false
    )

Get connected components of the given `system` based on a given cutoff for each
of the elements present, `nlcutoff` which by default is the covalent radii. 
`dists` and `retconnmat` specifies whether to return the distance matrix and 
connectivity matrix respectively.
"""
function getconnectedcomponents(
        str, dists=true;
        nlcutoff = natural_cutoffs(str).+0.2,
        retconnmat = false
    )   
    connmat, distmat = getconnectivitymatrix(
        str; nlcutoff=nlcutoff, retdistmat=true
    )
    conncomp = connectedcomponents(connmat)
    
    if dists
        if retconnmat
            return conncomp, distmat, connmat
        else
            return conncomp, distmat
        end
    else
        if retconnmat
            return conncomp, connmat
        else
            return conncomp
        end
    end
end

function connectedcomponents(connmat::AbstractMatrix)
    return connected_components(SimpleGraph(connmat))
end

"""
    transformpositions(f::Function, system::AbstractSystem)

Transform positions of the given `system` with `f` and return the transformed
system.
"""
# function transformpositions(f::Function, system::AbstractSystem)
#     positions = map(f, position(system))
#     cell = bounding_box(system)
#     at = species(system)
#     bc = boundary_conditions(system)
#     
#     typeof(system)(cell, positions, at)
# end
# 
function transformpositions(f::Function, system::FastSystem)
    particles = collect(system)
    positions = map(f, position.(particles))
    box = bounding_box(system)
    bc = boundary_conditions(system)
    
    FastSystem(box, bc, positions, atomic_symbol.(particles),
               atomic_number.(particles), atomic_mass.(particles))
end

"""
    wrap(system::AbstractSystem)

Wrap all the atoms in the given `system` into the cell and return the wrapped
system. Only works for orthogonal boxes.
"""
function wrap(system::AbstractSystem)
    origin = [0.0, 0.0, 0.0]u"Å"
    cell = box_lengths(system)
    f = posvec -> (ifelse.(
                           (origin .<= posvec .<= cell),
                           posvec,
                           posvec .+ (cell .* floor.((cell .- posvec) ./ cell))
                          )
    )

    return transformpositions(f, system)
end

"""
    interpolate_systems(sys1::AbstractSystem, sys2::AbstractSystem, nimg::Int)

Do a linear interpolation to obtain `nimg` structures between `sys1` and `sys2`.
"""
function interpolate_systems(
        sys1::FastSystem, sys2::FastSystem, nimg::Int
    )
    particles1 = collect(sys1)
    particles2 = collect(sys2)
    pos1 = position(sys1)
    pos2 = position(sys2)
    bcs = boundary_conditions(sys1)
    box = bounding_box(sys1)
    syms, nums, masses = atomic_symbol.(particles1), atomic_number.(particles1), atomic_mass.(particles1)

    if !(box == bounding_box(sys2))
        error("Can't interpolate systems with different cells.")
    end

    newsystems = Vector{typeof(sys1)}(undef, nimg)
    posstep = (pos2 .- pos1) / (nimg + 1)
    for i in 1:nimg
        newpos = pos1 .+ (i.*posstep)
        
        fsystem = FastSystem(box, bcs, newpos, syms, nums, masses)
        newsystems[i] = fsystem
    end

    return newsystems
end

end
