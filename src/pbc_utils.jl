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

function pbc_shortest_vectors(lattice::AbstractArray, fcoords1::AbstractArray, fcoords2::AbstractArray, return_dists=false, return_vects=true)
    I = size(fcoords1)[1]
    J = size(fcoords2)[1]
    T = eltype(lattice)

    # Allocates
    # cart_f1 = (fcoords1 .% 1) * lattice
    # cart_f2 = (fcoords2 .% 1) * lattice
    cart_f1 = dot_2d_mod(fcoords1, lattice)
    cart_f2 = dot_2d_mod(fcoords2, lattice)
    cart_im = images_view * lattice

    if return_vects
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
            if return_vects
                for l in 1:3
                    vectors[i, j, l] = pre_im[l] + cart_im[bestk, l]
                    vectors[j, i, l] = -vectors[i, j, l]
                end
            end
        end
    end

    if return_dists && return_vects
        return vectors, dists
    elseif return_dists
        return dists
    else
        return vectors
    end
end

function pbc_shortest_vectors(lattice::AbstractArray, fcoords::AbstractArray, return_dists=false, return_vects=true)
    I = size(fcoords)[1]
    T = eltype(lattice)

    cart_f = dot_2d_mod(fcoords, lattice)
    cart_im = images_view * lattice

    if return_vects
        vectors = similar(lattice, I, I, 3)
    end
    dists = similar(lattice, I, I)

    pre_im = similar(lattice, 3)
    for i in 1:I
        for j in i:I
            if i == j
                dists[i, j] = zero(T)
                if return_vects
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
            if return_vects
                for l in 1:3
                    vectors[i, j, l] = pre_im[l] + cart_im[bestk, l]
                    vectors[j, i, l] = -vectors[i, j, l]
                end
            end
        end
    end

    if return_dists && return_vects
        return vectors, dists
    elseif return_dists
        return dists
    else
        return vectors
    end
end
