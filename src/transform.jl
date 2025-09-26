function _getdefaultdata(system::AbstractSystem)
    particles = collect(system)
    positions = position(system)
    box = cell_vectors(system)
    bc = periodicity(system)

    return particles, positions, box, bc
end

"""
    transformpositions(f::Function, system::AbstractSystem)

Transform positions of the given `system` with `f` and return the transformed
system.
"""
# function transformpositions(f::Function, system::AbstractSystem)
#     positions = map(f, position(system))
#     cell = cell_vectors(system)
#     at = species(system)
#     bc = periodicity(system)
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
"""
function wrap(system)
    # 1) fractional (scaled) coordinates
    fracs = scaled_position(system)

    # 2) wrap each fractional triplet componentwise
    wrapped_fracs = [(f .- floor.(f)) for f in fracs]

    # 3) back to Cartesian
    newpos = cartesian_position.(Ref(system), wrapped_fracs)
    particles = Atom.(atomic_symbol.(system), newpos)

    FlexibleSystem(system; particles)
end
function wrap_position(sys, pos)
    f = scaled_position(sys, pos)
    wrapped_f = f .- floor.(f)
    wrapped_pos = cartesian_position(sys, wrapped_f)

    return wrapped_pos
end

"""
    supercell(system::AbstractSystem, supercellvec::Vector{Int}; sorted=false)

Create a supercell of the given `system`, where `supercellvec` is the
repetitions in each direction. All system and atom properties are copied to
the new system. If `sorted` is true, resultant supercell is sorted based on
atomic_symbol of input structure.
"""
function supercell(system::AbstractSystem, supercellvec::Vector{Int}; sorted=false)
    # Converts to fractional coordinates and adds unity in each of the direction
    # to be repeated
    if length(supercellvec) != 3 || !all(supercellvec .> 0)
        error("`supercellvec` should be a 3 element positive Vector.")
    end

    system_props = Dict{Symbol, Any}()
    # All system properties except cell_vectors is copied to the new system
    for (k,v) in pairs(system)
        if k != :cell_vectors
            system_props[k] = v
        end
    end
    
    cellmat = cell_matrix(system)
    newcellmat = cellmat .* supercellvec'
    system_props[:cell_vectors] = collect(eachcol(newcellmat))

    particles = collect(system)
    newparticles = copy(particles)

    # TODO sort based on atom type in input structure
    # newparticles = eltype(particles)[]

    # map(particles) do at
    #     spos = scaled_position(at, cellmat)

    #     map(1:3) do ax
    #         repunit = supercellvec[ax]
    #         for i in 0:(repunit-1)
    #             v = collect(spos)
    #             v[ax] += i
    #             
    #             at_props = Dict{Symbol, Any}()
    #             # All atom properties apart from position is copied
    #             for (k,v) in pairs(at)
    #                 if k != :position
    #                     at_props[k] = v
    #                 end
    #             end
    #             newpos = cellmat*v

    #             push!(newparticles, Atom(;position=newpos, at_props...))
    #         end
    #     end
    # end


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

    if sorted
        inp_atsyms = atomic_symbol(system)
        sort!(newparticles; by=x -> findfirst(==(atomic_symbol(x)), inp_atsyms))
    end

    FlexibleSystem(newparticles; system_props...)
end

"""
    sort(system::AbstractSystem; by=x->atomic_number(x), rev=false)

Sort an the atoms in an AbstractSystem. The `by` keyword argument specifies the
sorting function and takes as input the atom object. By default, the system is
sorted based on its atomic numbers. `rev` can be used to reverse the order.
"""
function sort(system::AbstractSystem; by=x->atomic_number(x), rev=false)
    atoms = collect(system)
    sortedatoms = sort(atoms; by, rev)
    return FlexibleSystem(sortedatoms; pairs(system)...)
end

"""
    reorder_new_system_by_reference(ref_sys, new_sys; tol=0.5u"Å")

- `ref_sys`: the original "reference" system (N atoms).
- `new_sys`: a modified system with M (MA) atoms added, so total N+M atoms.

Returns:
- A new PeriodicSystem with the **first N atoms** in the **same order** as `ref_sys`
  (matched by element and position, within `tol` distance),
- Followed by any unmatched atoms from `new_sys` (i.e., the newly added ones).

This ensures a 1-to-1 index correspondence for the old atoms between `ref_sys` and
the reordered `new_sys`.
"""
function reorder_new_system_by_reference(
    ref_sys::AbstractSystem,
    new_sys::AbstractSystem;
    tol = 0.5u"Å"
)
    # Extract relevant data for reference system
    ref_atoms = collect(ref_sys)
    ref_positions = position.(ref_atoms)
    ref_symbols   = atomic_symbol.(ref_atoms)

    # Extract data for new system
    new_atoms = collect(new_sys)
    new_positions = position.(new_atoms)
    new_symbols   = atomic_symbol.(new_atoms)

    # Keep track of atoms which are used up
    used = falses(length(new_atoms))

    # Helper function to find the closest new-atom index for a given ref-atom
    function find_closest_match(ref_sym, ref_pos; tol=tol)
        best_idx = nothing
        best_dist = Inf
        for (j, usedj) in enumerate(used)
            if !usedj && (new_symbols[j] == ref_sym)
                # same element => candidate
                # d = norm(ustrip.(u"Å", new_positions[j]) .- ustrip.(u"Å", ref_pos))
                # d = ustrip(u"Å", distance(ref_sys, new_positions[j], ref_pos))
                d = ustrip(u"Å", norm(AtomsToolbox.pbc_shortest_vector(ref_sys, new_positions[j], ref_pos)))
                if d < best_dist
                    best_dist = d
                    best_idx = j
                end
            end
        end
        if best_idx !== nothing && (isnothing(tol) || best_dist < ustrip(u"Å", tol))
            return best_idx
        else
            return nothing  # no suitable match found
        end
    end

    # Build array of new-system indices in the correct order
    matched_indices = Vector{Int}(undef, length(ref_atoms))

    # For each reference atom in order, find a matching new-atom
    for i in eachindex(ref_atoms)
        idx = find_closest_match(ref_symbols[i], ref_positions[i])
        if idx === nothing
            error("Could not find a match in new_sys for reference atom $i ($(ref_symbols[i])) within tolerance $tol")
        end
        matched_indices[i] = idx
        used[idx] = true
    end

    # Any new-atom not used in matching must be a newly added MA
    unmatched_indices = findall(!=(true), used)

    if !isempty(unmatched_indices) && !append_if_unmatched
        error("No match found for atoms: $(unmatched_indices). Set `append_if_unmatched` to append the atoms to the end.")
    else
        final_indices = vcat(matched_indices, unmatched_indices)
    end

    # List of Atoms in order
    reordered_atoms = Atom[]
    for idx in final_indices
        push!(reordered_atoms, new_atoms[idx])
    end

    return periodic_system(
        reordered_atoms,
        cell_vectors(new_sys);
        periodicity=periodicity(new_sys)
    )
end
