module CMInject

include("Particles.jl")
include("Fields.jl")

# ------ Implementation

function f!(du, u, p, t)
    ptype, pparams, fields = p
    particle = ptype(pparams, u)

    # TODO necessary to set to zero here? Or does DiffEq handle it for us?
    du[:] .= 0
    for field in fields
        du .+= acceleration(particle, field, t)
    end
    nothing
end

end