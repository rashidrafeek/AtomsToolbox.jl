# Functions to find shortest distance accounting for PBC, from pymatgen.
include("pbc_utils.jl")

"""
    cell_lengths(sys::AbstractSystem)

Obtain the cell lengths as a Vector, [a,b,c].
"""
cell_lengths(sys::AbstractSystem) = norm.(cell_vectors(sys))

"""
    cell_parameters(sys::AbstractSystem)

Obtain the cell lengths and angles as a Vector, [a,b,c,α,β,γ].
"""
function cell_parameters(sys::AbstractSystem)
    av,bv,cv = cell_vectors(sys)
    a,b,c = norm.((av,bv,cv))
    
    α = acosd((bv⋅cv)/(b*c))u"°"
    β = acosd((av⋅cv)/(a*c))u"°"
    γ = acosd((av⋅bv)/(a*b))u"°"

    return [a,b,c,α,β,γ]
end

"""
    cell_angles(sys::AbstractSystem)

Obtain the cell angles as a Vector, [α,β,γ].
"""
function cell_angles(sys::AbstractSystem)
    av,bv,cv = cell_vectors(sys)
    a,b,c = norm.((av,bv,cv))
    
    α = acosd((bv⋅cv)/(b*c))u"°"
    β = acosd((av⋅cv)/(a*c))u"°"
    γ = acosd((av⋅bv)/(a*b))u"°"

    return [α,β,γ]
end

"""
    distance(system::AbstractSystem, at1, at2)

Get the distance between atoms with indices `at1` and `at2`.
"""
function distance(system::AbstractSystem, at1::Int, at2::Int; pbc=all(periodicity(system)))
    pos1 = position(system, at1)
    pos2 = position(system, at2)

    return distance(system, pos1, pos2; pbc)
end
function distance(
        system::AbstractSystem, pos1::T, pos2::U; pbc=all(periodicity(system))
    ) where {T <: AbstractVector{<:Unitful.Length}, U <: AbstractVector{<:Unitful.Length}}
    if pbc
        dists = pbc_shortest_vectors(system, pos1, pos2, Val(true), Val(false))
        dist = only(dists)
    else
        dist = euclidean(pos1, pos2)
    end

    return dist
end

"""
    distance_matrix(system::AbstractSystem)

Get the distance matrix between all atoms as an NxN matrix where N is the
number of atoms in the given `system`.
"""
function distance_matrix(system::AbstractSystem; pbc=all(periodicity(system)))
    if pbc
        dists = pbc_shortest_vectors(system, position(system, :), Val(true), Val(false))
    else
        dists = pairwise(Euclidean(), position(system, :))
    end

    return dists
end

"""
    angle(system::AbstractSystem, at1, at2, at3)

Get the angle between vectors connecting atoms with indices `at2`-`at1` and 
`at2`-`at3`.
"""
function angle(system::AbstractSystem, at1::Int, at2::Int, at3::Int; pbc=all(periodicity(system)))
    pos1 = position(system, at1)
    pos2 = position(system, at2)
    pos3 = position(system, at3)

    return angle(system, pos1, pos2, pos3; pbc)
end
function angle(
        system::AbstractSystem, pos1::T, pos2::T, pos3::T; pbc=all(periodicity(system))
    ) where {T <: AbstractVector{<:Unitful.Length}}
    if pbc
        cell = reduce(hcat, cell_vectors(system))'
        icell = inv(cell)
        frpos1 = pos1' * icell
        frpos2 = pos2' * icell
        frpos3 = pos3' * icell
        vec1 = pbc_shortest_vectors(cell, frpos2, frpos1)
        vec2 = pbc_shortest_vectors(cell, frpos2, frpos3)
    else
        vec1 = pos2 - pos1
        vec2 = pos2 - pos3
    end

    cos_ang = (vec1 ⋅ vec2)/(norm(vec1)*norm(vec2))
    ang = abs(cos_ang) <= one(cos_ang) ? acosd(cos_ang)u"°" : _slow_acosd(cos_ang)u"°"
    # ang = acosd(cos_ang)u"°"

    return ang
end
# For edge case when x >= 1. See https://discourse.julialang.org/t/linearalgebra-unstable-floating-point-operations/45581/5
@noinline _slow_acosd(x) = x ≈ 1 ? zero(x) : x ≈ -1 ? one(x)*180 : acosd(x)

"""
    dihedral(system::AbstractSystem, at1, at2, at3, at4)

Get the dihedral angle for the atoms specified by indices `at1`, `at2`, `at3`
and `at4`.
"""
function dihedral(system::AbstractSystem, at1::Int, at2::Int, at3::Int, at4::Int; pbc=all(periodicity(system)))
    pos1 = position(system, at1)
    pos2 = position(system, at2)
    pos3 = position(system, at3)
    pos4 = position(system, at4)

    return dihedral(system, pos1, pos2, pos3, pos4; pbc)
end
function dihedral(
        system::AbstractSystem, pos1::T, pos2::T, pos3::T, pos4::T; pbc=all(periodicity(system))
    ) where {T <: AbstractVector{<:Unitful.Length}}
    # Source: Blondel and Karplus, J. Comp. Chem., Vol. 17, No. 9, 1
    #   132-1 141 (1 996).
    # if pbc
    #     cell = reduce(hcat, cell_vectors(system))'
    #     icell = inv(cell)
    #     frpos1 = pos1' * icell
    #     frpos2 = pos2' * icell
    #     frpos3 = pos3' * icell
    #     frpos4 = pos4' * icell
    #     vec0 = -pbc_shortest_vectors(cell, frpos3, frpos2)[1,1,:]  # Axis
    #     vec1 = pbc_shortest_vectors(cell, frpos2, frpos1)[1,1,:]
    #     vec2 = pbc_shortest_vectors(cell, frpos4, frpos3)[1,1,:]
    # else
    #     vec0 = pos2 - pos3
    #     vec1 = pos2 - pos1
    #     vec2 = pos4 - pos3
    # end
    # p1 = vec1 × vec0
    # p2 = vec0 × vec2
    # cross_p1p2 = p1 × p2

    # sin_ang = (cross_p1p2 ⋅ vec0)/(norm(p1)*norm(p2)*norm(vec0))
    # cos_ang = (p1 ⋅ p2)/(norm(p1)*norm(p2))

    # ang = atand(sin_ang, cos_ang)u"°"
    
    # Obtained from pymatgen
    if pbc
        cell = reduce(hcat, cell_vectors(system))'
        icell = inv(cell)
        frpos1 = pos1' * icell
        frpos2 = pos2' * icell
        frpos3 = pos3' * icell
        frpos4 = pos4' * icell
        vec1 = pbc_shortest_vectors(cell, frpos1, frpos2)[1,1,:]  # Axis
        vec2 = pbc_shortest_vectors(cell, frpos2, frpos3)[1,1,:]
        vec3 = pbc_shortest_vectors(cell, frpos3, frpos4)[1,1,:]
    else
        vec1 = pos2 - pos1
        vec2 = pos3 - pos2
        vec3 = pos4 - pos3
    end
    p23 = vec2 × vec3
    p12 = vec1 × vec2

    ang = atand(norm(vec2)*(vec1 ⋅ p23), p12 ⋅ p23)u"°"

    # if ang < 0.0u"°"
    #     return 360u"°" + ang
    # else
    #     return ang
    # end
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
    connectivity_matrix(
        system;                                                                 
        nlcutoff = natural_cutoffs(system).+0.2,
        retdistmat = false
    )

Get the connectivity matrix of the given `system`. The cutoff is specified
using `nlcutoff`. To return the distance matrix, set `retdistmat` as `true`.
"""
function connectivity_matrix(system::AbstractSystem;
                               nlcutoff=natural_cutoffs(system) .+ 0.2u"Å",
                               retdistmat=false)
    distmat = distance_matrix(system)

    connmat = (distmat .- nlcutoff') .< nlcutoff

    if retdistmat
        return connmat, distmat
    else
        return connmat
    end
end
function connectivity_matrix(distmat::Matrix; nlcutoff)
    connmat = (distmat .- nlcutoff') .< nlcutoff

    return connmat
end

"""
    connected_components(
        system, dists=true;
        nlcutoff = natural_cutoffs(str).+0.2,
        retconnmat = false
    )

Get connected components of the given `system` based on a given cutoff for each
of the elements present, `nlcutoff` which by default is the covalent radii. 
`dists` and `retconnmat` specifies whether to return the distance matrix and 
connectivity matrix respectively.
"""
function connected_components(
        str,
        dists=true;
        nlcutoff=natural_cutoffs(str) .+ 0.2u"Å",
        retconnmat=false
    )
    connmat, distmat = connectivity_matrix(str; nlcutoff=nlcutoff, retdistmat=true)
    conncomp = connected_components(connmat)

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

function connected_components(connmat::AbstractMatrix)
    return Graphs.connected_components(Graphs.SimpleGraph(connmat))
end

"""
    scaled_position(sys::AbstractSystem)
    scaled_position(sys::AbstractSystem, index)
    scaled_position(sys::AbstractSystem, pos::AbstractVector)
    scaled_position(atom::Union{Atom,AtomView}, cellmat)
    scaled_position(pos::AbstractVector{<:Unitful.Length}, cellmat)

Obtain the scaled positions (fractional coordinates) with respect to the cell
matrix for the given system or `atom`.
"""
function scaled_position(sys::AbstractSystem)
    cellmat = cell_matrix(sys)
    frpos = reduce(hcat, position(sys, :))' * inv(cellmat')

    Vector.(eachrow(frpos))
end
scaled_position(sys::AbstractSystem, index) = scaled_position(sys[index], cell_matrix(sys))
scaled_position(sys::AbstractSystem, pos::AbstractVector) = scaled_position(pos, cell_matrix(sys))
scaled_position(atom::Union{Atom,AtomView}, cellmat) = scaled_position(position(atom), cellmat)
function scaled_position(pos::AbstractVector{<:Unitful.Length}, cellmat)
    frpos = pos' * inv(cellmat')

    frpos'
end

"""
    cartesian_position(sys::AbstractSystem, pos::AbstractVector)
    cartesian_position(atom::Atom, cellmat)
    cartesian_position(pos::AbstractVector, cellmat)

Obtain the Cartesian positions with respect to the cell matrix for the given system, `atom`, or `pos`.
"""
cartesian_position(sys::AbstractSystem, pos::AbstractVector) = cartesian_position(pos, cell_matrix(sys))
cartesian_position(atom::Union{Atom, AtomView}, cellmat) = cartesian_position(position(atom), cellmat)
function cartesian_position(pos::AbstractVector{<:Real}, cellmat)
    cart_pos = pos' * cellmat'  # Convert fractional positions to Cartesian coordinates
    cart_pos'
end

#
# Cell utils
# 

@inline cell_matrix(sys::AbstractSystem) = reduce(hcat, cell_vectors(sys))
# lu is required as unitful matrices are erroring without it
# when the cell matrix is diagonal
"""
    volume(sys::AbstractSystem)

Obtain the volume of the given system.
"""
volume(sys::AbstractSystem) = det(lu(cell_matrix(sys)))

# 
# Generic utils
#

# PBC for non-rectangular cells
# See: https://github.com/JuliaStats/Distances.jl/pull/237
# function Distances._evaluate(
#         dist::PeriodicEuclidean, ua::AbstractVector{Q}, ub::AbstractVector{Q}, up::AbstractMatrix{Q}
#     ) where {Q <: Unitful.AbstractQuantity}
#     a,b,p = ustrip(ua), ustrip(ub), ustrip(up) # Currently solve doesn't work with Unitful matrices
#     ndim = LinearAlgebra.checksquare(p)
#     @boundscheck if length(a) != length(b)
#         throw(DimensionMismatch("first array has length $(length(a)) which does not match the length of the second, $(length(b))."))
#     end
#     @boundscheck if length(a) != ndim
#         throw(DimensionMismatch("arrays have length $(length(a)) but basis vectors have length $ndim."))
#     end
#     s = p \ (a - b)      # decompose to p-basis, a - b = sx x⃗ + sy y⃗ + ...
#     mindistance = Inf
#     s = mod.(s, 1)       # move to first "quadrant"
#     @inbounds for k = 0:1<<ndim-1  # 2^{# of dimensions} possible choices
#         loc = sum(i->view(p,:,i) .* (k>>(i-1) & 1 == 0 ? s[i] : s[i]-1), 1:ndim)
#         mindistance = min(norm(loc), mindistance)
#     end
#     return mindistance*unit(Q)
# end
