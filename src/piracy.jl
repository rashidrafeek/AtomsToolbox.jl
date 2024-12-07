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

# These were removed in the AtomsBase 0.4 update
Base.position(sys::AbstractSystem) = position.(sys, :)
AtomsBase.atomic_symbol(sys::AbstractSystem) = atomic_symbol.(species(sys, :))
AtomsBase.atomic_number(sys::AbstractSystem) = atomic_number.(species(sys, :))
