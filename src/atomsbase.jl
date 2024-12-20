# """
#     _getconstructor(s::AbstractSystem)
# 
# Get the constructor to construct the concrete type of `s` which is a subtype of
# AbstractSystem
# """
# function _getconstructor(s::AbstractSystem)
#     return error("""
# No constructor defined. Please define _getconstructor(::$(Base.typename(typeof(s)).name)).
# """)
# end
# _getconstructor(::FastSystem) = FastSystem
# _getconstructor(::FlexibleSystem) = FlexibleSystem
# 
# function Base.getindex(system::AbstractSystem, range::AbstractVector)
#     particles = collect(system)
#     return _getconstructor(system)(particles[range],
#                                    cell_vectors(system),
#                                    periodicity(system))
#     # _getconstructor(system)(
#     #     cell_vectors(system), periodicity(system), 
#     #     position(system)[range], atomic_symbol(system)[range],
#     #     atomic_number(system)[range], atomic_mass(system)[range])
# end
