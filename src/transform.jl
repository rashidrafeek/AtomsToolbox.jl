function _getdefaultdata(system::AbstractSystem)
    particles = collect(system)
    positions = position(system)
    box = bounding_box(system)
    bc = boundary_conditions(system)

    return particles, positions, box, bc
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
    particles, positions, box, bc = _getdefaultdata(system)
    final_positions = map(f, positions)

    return FastSystem(box,
                      bc,
                      final_positions,
                      atomic_symbol.(particles),
                      atomic_number.(particles),
                      atomic_mass.(particles))
end

"""
    wrap(system::AbstractSystem)

Wrap all the atoms in the given `system` into the cell and return the wrapped
system. 

!!! warning
    Currently only works for orthogonal boxes.
"""
function wrap(system::AbstractSystem)
    origin = [0.0, 0.0, 0.0]u"Å"
    cell = cell_lengths(system)
    f = posvec -> (ifelse.((origin .<= posvec .<= cell),
                           posvec,
                           posvec .+ (cell .* floor.((cell .- posvec) ./ cell))))

    # cellmat = getcellmatrix(system)'
    # icell = inv(cellmat)
    # f = function (posvec)
    #     frpos = posvec * icell
    #     frposincell = sign.(frpos) .* getindex.(modf.(frpos),1)
    #     return frposincell*cellmat'
    # end
    return transformpositions(f, system)
end

# function makesupercell(system::FastSystem, supercellmatrix)
#     particles, positions, box, bc = _getdefaultdata(system)
# 
#     
# 
#     return FastSystem(box,
#                       bc,
#                       positions,
#                       atomic_symbol.(particles),
#                       atomic_number.(particles),
#                       atomic_mass.(particles))
# end
