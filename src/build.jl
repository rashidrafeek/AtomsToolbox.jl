"""
    interpolate_systems(sys1::AbstractSystem, sys2::AbstractSystem, nimg::Int)

Do a linear interpolation to obtain `nimg` structures between `sys1` and `sys2`.
"""
function interpolate_systems(sys1::FastSystem, sys2::FastSystem, nimg::Int)
    particles1 = collect(sys1)
    particles2 = collect(sys2)
    pos1 = position(sys1)
    pos2 = position(sys2)
    bcs = periodicity(sys1)
    box = cell_vectors(sys1)
    syms, nums, masses = atomic_symbol.(particles1),
                         atomic_number.(particles1),
                         atomic_mass.(particles1)

    if !(box == cell_vectors(sys2))
        error("Can't interpolate systems with different cells.")
    end

    newsystems = Vector{typeof(sys1)}(undef, nimg)
    posstep = (pos2 .- pos1) / (nimg + 1)
    for i in 1:nimg
        newpos = pos1 .+ (i .* posstep)

        fsystem = FastSystem(box, bcs, newpos, syms, nums, masses)
        newsystems[i] = fsystem
    end

    return newsystems
end
