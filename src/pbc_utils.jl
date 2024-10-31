# Translation of coord_cython.pyx from pymatgen
# See https://github.com/materialsproject/pymatgen/blob/master/pymatgen/util/coord_cython.pyx

const images_view = permutedims(reduce(hcat, 
                   [[a,b,c] for (a,b,c) in Iterators.product(-1:1, -1:1, -1:1)][:])
               )

function dot_2d_mod(a,b)
    I = size(a)[1]
    J = size(b)[2]
    K = size(a)[2]
    T = eltype(b)
    
    out = zeros(T, I, J)
    for j in 1:J
        for k in 1:K
            for i in 1:I
                out[i,j] += a[i,k] % 1 * b[k,j]
            end
        end
    end
    
    return out
end

"""
    pbc_shortest_vectors(
        system::AbstractSystem, at1::Int, at2::Int,
        return_dists::Val{RD}=Val(false), return_vects::Val{RV}=Val(true)
    ) where {RD, RV}
    pbc_shortest_vectors(
        system::AbstractSystem, pos1::T, pos2::T, 
        return_dists::Val{RD}=Val(false), return_vects::Val{RV}=Val(true)
    ) where {T <: AbstractVector{<: Unitful.Length}, RD, RV}

Obtain the shortest vectors between position vectors, `pos1` and `pos2`,
taking into account PBC of the given `system`.
"""
function pbc_shortest_vectors(
        system::AbstractSystem, at1::Int, at2::Int,
        return_dists::Val{RD}=Val(false), return_vects::Val{RV}=Val(true)
    ) where {RD, RV}
    pos1 = position(system, at1)
    pos2 = position(system, at2)

    return pbc_shortest_vectors(system, pos1, pos2, return_dists, return_vects)[1,1,:]
end
function pbc_shortest_vectors(
        system::AbstractSystem, pos1::T, pos2::U,
        return_dists::Val{RD}=Val(false), return_vects::Val{RV}=Val(true)
    ) where {T <: AbstractVector{<: Unitful.Length}, U <: AbstractVector{<: Unitful.Length}, RD, RV}
    cell = reduce(hcat, bounding_box(system))'
    icell = inv(cell)
    frpos1 = pos1' * icell
    frpos2 = pos2' * icell

    return pbc_shortest_vectors(cell, frpos1, frpos2, return_dists, return_vects)
end
function pbc_shortest_vectors(
        system::AbstractSystem, posvec::T,
        return_dists::Val{RD}=Val(false), 
        return_vects::Val{RV}=Val(true)
    ) where {T <: Vector{<: AbstractVector{<: Unitful.Length}}, RD, RV}
    cell = reduce(hcat, bounding_box(system))'
    pos = reduce(hcat, posvec)'
    frpos = pos * inv(cell)

    return pbc_shortest_vectors(cell, frpos, return_dists, return_vects)
end
function pbc_shortest_vectors(
        lattice::AbstractArray, 
        fcoords1::AbstractArray, fcoords2::AbstractArray,
        return_dists::Val{RD}=Val(false), 
        return_vects::Val{RV}=Val(true)
    ) where {RD,RV}
    I = size(fcoords1)[1]
    J = size(fcoords2)[1]
    T = eltype(lattice)

    # Allocates
    # cart_f1 = (fcoords1 .% 1) * lattice
    # cart_f2 = (fcoords2 .% 1) * lattice
    cart_f1 = dot_2d_mod(fcoords1, lattice)
    cart_f2 = dot_2d_mod(fcoords2, lattice)
    cart_im = images_view * lattice

    if RV
        vectors = similar(lattice, I, J, 3)
    end
    dists = similar(lattice, I, J)

    pre_im = similar(lattice, 3)
    for i in 1:I
        for j in 1:J
            for l in 1:3
                pre_im[l] = cart_f2[j,l] .- cart_f1[i,l]
            end
            best = 1e100*oneunit(T)^2
            bestk = 100
            for k in 1:27
                da = (pre_im[1] + cart_im[k, 1])^2
                db = (pre_im[2] + cart_im[k, 2])^2
                dc = (pre_im[3] + cart_im[k, 3])^2
                d = da + db + dc
                if d < best
                    best = d
                    bestk = k
                end
            end
            dists[i,j] = sqrt(best)
            if RV
                for l in 1:3
                    vectors[i, j, l] = pre_im[l] + cart_im[bestk, l]
                    # vectors[j, i, l] = -vectors[i, j, l]
                end
            end
        end
    end

    if RD && RV
        return vectors, dists
    elseif RD
        return dists
    else
        return vectors
    end
end

function pbc_shortest_vectors(
        lattice::AbstractArray, 
        fcoords::AbstractArray, 
        return_dists::Val{RD}=Val(false), 
        return_vects::Val{RV}=Val(true)
    ) where {RD,RV}
    I = size(fcoords)[1]
    T = eltype(lattice)

    cart_f = dot_2d_mod(fcoords, lattice)
    cart_im = images_view * lattice

    if RV
        vectors = similar(lattice, I, I, 3)
    end
    dists = similar(lattice, I, I)

    pre_im = similar(lattice, 3)
    for i in 1:I
        for j in i:I
            if i == j
                dists[i, j] = zero(T)
                if RV
                    vectors[i, j, :] .= zero(T)
                end
                continue
            end
            for l in 1:3
                pre_im[l] = cart_f[j,l] - cart_f[i,l]
            end
            best = 1e100*oneunit(T)^2
            bestk = 100
            for k in 1:27
                da = (pre_im[1] + cart_im[k, 1])^2
                db = (pre_im[2] + cart_im[k, 2])^2
                dc = (pre_im[3] + cart_im[k, 3])^2
                d = da + db + dc
                if d < best
                    best = d
                    bestk = k
                end
            end
            dists[i,j] = dists[j,i] = sqrt(best)
            if RV
                for l in 1:3
                    vectors[i, j, l] = pre_im[l] + cart_im[bestk, l]
                    vectors[j, i, l] = -vectors[i, j, l]
                end
            end
        end
    end

    if RD && RV
        return vectors, dists
    elseif RD
        return dists
    else
        return vectors
    end
end

function pbc_shortest_vectors_withpbcimg(
        lattice::AbstractArray, 
        fcoords::AbstractArray, 
        return_dists::Val{RD}=Val(false), 
        return_vects::Val{RV}=Val(true)
    ) where {RD,RV}
    I = size(fcoords)[1]
    T = eltype(lattice)

    cart_f = dot_2d_mod(fcoords, lattice)
    cart_im = images_view * lattice

    if RV
        vectors = similar(lattice, I, I, 3)
    end
    dists = similar(lattice, I, I)

    pre_im = similar(lattice, 3)
    pbcimages = Array{Tuple{Int,Int,Int}}(undef, I, I)
    for i in 1:I
        for j in i:I
            if i == j
                dists[i, j] = zero(T)
                pbcimages[i,j] = (0,0,0)
                if RV
                    vectors[i, j, :] .= zero(T)
                end
                continue
            end
            for l in 1:3
                pre_im[l] = cart_f[j,l] - cart_f[i,l]
            end
            best = 1e100*oneunit(T)^2
            bestk = 100
            for k in 1:27
                da = (pre_im[1] + cart_im[k, 1])^2
                db = (pre_im[2] + cart_im[k, 2])^2
                dc = (pre_im[3] + cart_im[k, 3])^2
                d = da + db + dc
                if d < best
                    best = d
                    bestk = k
                end
            end
            dists[i,j] = dists[j,i] = sqrt(best)
            pbcimages[i,j] = Tuple(images_view[bestk,:])
            pbcimages[j,i] = Tuple(-images_view[bestk,:])
            if RV
                for l in 1:3
                    vectors[i, j, l] = pre_im[l] + cart_im[bestk, l]
                    vectors[j, i, l] = -vectors[i, j, l]
                end
            end
        end
    end

    return vectors, dists, pbcimages
end
