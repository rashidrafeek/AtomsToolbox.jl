box_lengths(sys::AbstractSystem) = getindex.(bounding_box(sys), [1, 2, 3])

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
        pairwise(PeriodicEuclidean(ustrip.(box_lengths(system))), ustrip.(position(system)))
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
function getconnectivitymatrix(system;
                               nlcutoff=natural_cutoffs(system) .+ 0.2,
                               retdistmat=false)
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
function getconnectedcomponents(str,
                                dists=true;
                                nlcutoff=natural_cutoffs(str) .+ 0.2,
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
