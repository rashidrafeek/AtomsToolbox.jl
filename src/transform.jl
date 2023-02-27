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

"""
    supercell(system::AbstractSystem, supercellvec::Vector{Int})

Create a supercell of the given `system`, where `supercellvec` is the
repetitions in each direction. All system and atom properties are copied to
the new system.
"""
function supercell(system::AbstractSystem, supercellvec::Vector{Int})
    # Converts to fractional coordinates and adds unity in each of the direction
    # to be repeated
    if length(supercellvec) != 3 || !all(supercellvec .> 0)
        error("`supercellvec` should be a 3 element positive Vector.")
    end

    system_props = Dict{Symbol, Any}()
    # All system properties except bounding_box is copied to the new system
    for (k,v) in pairs(system)
        if k != :bounding_box
            system_props[k] = v
        end
    end
    
    cellmat = cell_matrix(system)
    newcellmat = cellmat .* supercellvec'
    system_props[:bounding_box] = collect(eachcol(newcellmat))

    particles = collect(system)
    newparticles = copy(particles)
    for ax in 1:3
        repunit = supercellvec[ax]

        append!(newparticles, map(1:(repunit-1)) do i
            map(newparticles) do at
                v = collect(scaled_position(at, cellmat))
                v[ax] += i
                
                at_props = Dict{Symbol, Any}()
                # All atom properties apart from position is copied
                for (k,v) in pairs(at)
                    if k != :position
                        at_props[k] = v
                    end
                end
                newpos = cellmat*v

                Atom(;position=newpos, at_props...)
            end
        end...)
    end

    FlexibleSystem(newparticles; system_props...)
end
