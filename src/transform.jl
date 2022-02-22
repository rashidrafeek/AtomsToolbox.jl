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

    return FastSystem(box,
                      bc,
                      positions,
                      atomic_symbol.(particles),
                      atomic_number.(particles),
                      atomic_mass.(particles))
end

"""
    wrap(system::AbstractSystem)

Wrap all the atoms in the given `system` into the cell and return the wrapped
system. Only works for orthogonal boxes.
"""
function wrap(system::AbstractSystem)
    origin = [0.0, 0.0, 0.0]u"Å"
    cell = box_lengths(system)
    f = posvec -> (ifelse.((origin .<= posvec .<= cell),
                           posvec,
                           posvec .+ (cell .* floor.((cell .- posvec) ./ cell))))

    return transformpositions(f, system)
end
