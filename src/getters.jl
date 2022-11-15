# Workaround for Unitful StaticMatrix not working.
# See: https://github.com/PainterQubits/Unitful.jl/issues/538
function Base.inv(x::StaticMatrix{N,M,T}) where {N,M,T <: Unitful.AbstractQuantity}
    m = inv(ustrip.(x))
    iq = eltype(m)
    reinterpret(
        Unitful.Quantity{iq, inv(Unitful.dimension(T)), typeof(inv(unit(T)))}, 
        m
    )
end 

# Functions to find shortest distance accounting for PBC, from pymatgen.
include("pbc_utils.jl")

"""
    cell_lengths(sys::AbstractSystem)

Obtain the cell lengths as a Vector, [a,b,c].
"""
cell_lengths(sys::AbstractSystem) = norm.(bounding_box(sys))

"""
    cell_lengths_and_angles(sys::AbstractSystem)

Obtain the cell lengths and angles as a Vector, [a,b,c,α,β,γ].
"""
function cell_lengths_and_angles(sys::AbstractSystem)
    av,bv,cv = bounding_box(sys)
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
    av,bv,cv = bounding_box(sys)
    a,b,c = norm.((av,bv,cv))
    
    α = acosd((bv⋅cv)/(b*c))u"°"
    β = acosd((av⋅cv)/(a*c))u"°"
    γ = acosd((av⋅bv)/(a*b))u"°"

    return [α,β,γ]
end

"""
    getdistance(system::AbstractSystem, at1, at2)

Get the distance between atoms with indices `at1` and `at2`.
"""
function getdistance(system::AbstractSystem, at1, at2)
    pos1 = position(system, at1)
    pos2 = position(system, at2)
    if all(periodicity(system))
        cell = reduce(hcat, bounding_box(system))'
        icell = inv(cell)
        frpos1 = pos1' * icell
        frpos2 = pos2' * icell
        dists = pbc_shortest_vectors(cell, frpos1, frpos2, Val(true), Val(false))
        dist = only(dists)
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
function getdistancematrix(system::AbstractSystem, pbc=nothing)
    if isnothing(pbc)
        if all(periodicity(system))
            cell = reduce(hcat, bounding_box(system))'
            pos = reduce(hcat, position(system))'
            frpos = pos * inv(cell)

            dists = pbc_shortest_vectors(cell, frpos, Val(true), Val(false))
        else
            dists = pairwise(Euclidean(), position(system))
        end

        return dists
    else
        if pbc
            cell = reduce(hcat, bounding_box(system))'
            pos = reduce(hcat, position(system))'
            frpos = pos * inv(cell)

            dists = pbc_shortest_vectors(cell, frpos, Val(true), Val(false))
        else
            dists = pairwise(Euclidean(), position(system))
        end

        return dists
    end
end

"""
    getangle(system::AbstractSystem, at1, at2, at3)

Get the angle between vectors connecting atoms with indices `at2`-`at1` and 
`at2`-`at3`.
"""
function getangle(system::AbstractSystem, at1, at2, at3)
    pos1 = position(system, at1)
    pos2 = position(system, at2)
    pos3 = position(system, at3)
    if all(periodicity(system))
        cell = reduce(hcat, bounding_box(system))'
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

    ang = acosd((vec1 ⋅ vec2)/(norm(vec1)*norm(vec2)))u"°"

    return ang
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
function getconnectivitymatrix(system::AbstractSystem;
                               nlcutoff=natural_cutoffs(system) .+ 0.2u"Å",
                               retdistmat=false)
    distmat = getdistancematrix(system)

    connmat = (distmat .- nlcutoff') .< nlcutoff

    if retdistmat
        return connmat, distmat
    else
        return connmat
    end
end
function getconnectivitymatrix(distmat::Matrix; nlcutoff)
    connmat = (distmat .- nlcutoff') .< nlcutoff

    return connmat
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
function getconnectedcomponents(str,
                                dists=true;
                                nlcutoff=natural_cutoffs(str) .+ 0.2u"Å",
                                retconnmat=false)
    connmat, distmat = getconnectivitymatrix(str; nlcutoff=nlcutoff, retdistmat=true)
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

function getfractionalcoordinates(sys::AbstractSystem)
    cellmat = getcellmatrix(sys)
    frpos = reduce(hcat, position(sys))' * inv(cellmat')
end

#
# Cell utils
# 

@inline getcellmatrix(sys::AbstractSystem) = reduce(hcat, bounding_box(sys))
getvolume(sys::AbstractSystem) = det(getcellmatrix(sys))

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
