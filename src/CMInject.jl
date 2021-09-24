module CMInject

include("Particles.jl")
include("Fields.jl")

using DifferentialEquations

# ------ Implementation

function f!(du, u, params, t)
    p, fields = params
    particle = SphericalParticle{Float64}(u, p)

    du .= 0
    du.x = u.vx
    du.z = u.vz
    for field in fields
        ax, az = acceleration(particle, field, t)
        du.vx = ax
        du.vz = az
    end
    nothing
end

function f_direct!(du, u, params, t)
    p, fields = params
    particle = SphericalParticle{Float64}(u, p)

    du .= 0
    du.x = u.vx
    du.z = u.vz
    for field in fields
        ax, az = acceleration_direct(particle, field, t)
        du.vx = ax
        du.vz = az
    end
    nothing
end

const du = similar(CMInject.example_particle._u)
function quick()
    p = CMInject.example_particle
    f = CMInject.example_field
    f!(du, p._u, (p._p, (f,)), 0.0)
    du
end

function quick_direct()
    p = CMInject.example_particle
    f = CMInject.example_field
    f_direct!(du, p._u, (p._p, (f,)), 0.0)
    du
end

@inline function always_false(args...)
    false
end

function problem(x0, z0)
    p = CMInject.example_particle
    fs = (CMInject.example_field,)
    u0 = similar(p._u)
    u0 .= p._u
    u0.x = x0
    u0.z = z0
    prob = ODEProblem{true}(f!, u0, (0.0, 0.03), (p._p, fs))
    sol = solve(prob, Heun(), dt=1e-5, adaptive=false, dense=false, unstable_check=always_false)
    sol
end

end