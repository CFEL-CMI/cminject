module CMInject

include("Particles.jl")
include("Fields.jl")

# ------ Implementation

function f!(du, u, p, t)
    pparams, field = p
    particle = SphericalParticle(pparams, u)

    # TODO necessary to set to zero here? Or does DiffEq handle it for us?
    # TODO make sure this only applies to the velocity components...
    @inbounds begin
        du .= 0
        du.x = du.vx
        du.z = du.vz

        ax, az = acceleration(particle, field, t)
        du.vx += ax
        du.vz += az
    end
    nothing
end

end