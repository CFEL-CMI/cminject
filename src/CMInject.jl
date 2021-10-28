module CMInject

include("Interpolation.jl")
include("Particles.jl")
include("Fields.jl")
include("Sources.jl")
include("Detectors.jl")
include("ResultPlots.jl")

using DifferentialEquations


# Generic problem/experiment implementation

@generated function f!(du, u, params, t)
    # TODO this isn't so nice, we need to pass an initial p twice:
    #  - for determining particle type (indirectly, via associated_particle_type on p's type at compile time)
    #  - the actual p that can be modified, during the integration loop
    #    -> we can't pass this on the type-level since it will live in a mutable list (I think)
    ptype = associated_particle_type(params.types[1])
    quote
        _, (p,), fields = params
        particle = $ptype(u, p)

        du .= 0
        # TODO refactor to be generic, based on the @velocities meta-information from a future @particle macro
        transfer_velocities!(du, particle)
        @inbounds for field in fields
            # TODO refactor to be generic, based on the values in the @velocities meta-information
            acc = acceleration(particle, field, time)
            accelerate!(du, acc)
        end
        nothing
    end
end

# TODO the same todo notes as for f! apply here
@generated function g!(du, u, params, t)
    ptype = associated_particle_type(params.types[1])
    quote
        _, (p,), fields = params
        particle = $ptype(u, p)

        du .= 0
        @inbounds for field in fields
            acc = noise(particle, field, t)
            accelerate!(du, acc)
        end
        nothing
    end
end


# Initial values helper fn for the toy ODE/SDE problems below

function initial_vals(x0, z0)
    p = CMInject.example_particle
    u0 = similar(p._u)
    u0 .= p._u
    u0.x = x0
    u0.z = z0
    fs = (CMInject.example_field,)
    params = (p._p, [p._p], fs)
    tspan = (0.0, 0.03)
    u0, tspan, params
end


# Default options for the solvers

@inline function always_false(args...)
    false
end

const default_solver_opts = (dt=1e-5, adaptive=false, dense=false, unstable_check=always_false)


## ODE problem

function ode_problem(x0, z0)
    u0, tspan, params = initial_vals(x0, z0)
    ODEProblem{true}(f!, u0, tspan, params)
end

function ode_solve(args...)
    prob = ode_problem(args...)
    solver = Heun()
    solve(prob, solver; default_solver_opts...)
end


## SDE problem

function sde_problem(x0, z0)
    u0, tspan, params = initial_vals(x0, z0)
    SDEProblem{true}(f!, g!, u0, tspan, params)
end

function sde_solve(args...)
    prob = sde_problem(args...)
    solver = EulerHeun()
    solve(prob, solver; default_solver_opts...)
end


## SDE Ensemble problem

const example_dists = Dict(
    :x => Normal(0, 1e-3), :z => Dirac(-0.128),
    :vx => Normal(0, 0.1), :vz => Dirac(10.0),
    :r => Normal(13.5e-9, 2e-9), :ρ => Dirac(1050.0)
)
const example_source = SamplingSource{SphericalParticle2D{Float64}}(example_dists)

function make_initial_values(source, fields, n)
    particles = generate(source, n)
    function prob_func(prob::P, i, repeat)::P where P
        particle = particles[i]
        prob.u0 .= particle._u
        prob.p[2][1] = particle._p
        prob
    end

    u0 = similar(particles[1]._u)
    u0 .= particles[1]._u
    p0 = particles[1]._p
    params = (p0, [p0], fields)

    u0, params, prob_func
end

function simulate(
    ; source::S, fields, detectors, tspan, n_particles::Integer,
    solver=EulerHeun(), solver_opts=(), ensemble_strategy=EnsembleThreads()
) where {S<:Source}
    make_initial_values(source, fields, n_particles)
    u0, params, prob_func = make_initial_values(source, fields, n_particles)
    problem = SDEProblem{true}(f!, g!, u0, tspan, params)
    ep = EnsembleProblem(problem, prob_func=prob_func, safetycopy=false)
    sol = solve(ep, solver, ensemble_strategy; trajectories=n_particles, solver_opts...)

    detector_hits = [calculate_hits(sol, detector) for detector in detectors]
    return sol, detector_hits
end


function ensemble_solve(n; tspan=(0, 0.05), detector_positions=-0.003:0.001:0.003)
    simulate(
        ; source=example_source, fields=(example_field,),
        detectors=Tuple([SectionDetector{Float64, :z}(pos, true) for pos in detector_positions]),
        tspan=tspan, n_particles=n,
        solver_opts=default_solver_opts
    )
end


## Benchmarking utilities - not to actually solve physical problems

const _du = similar(CMInject.example_particle._u)
function _quick()
    p = CMInject.example_particle
    f = CMInject.example_field
    f!(_du, p._u, (p._p, [p._p,], (f,)), 0.0)
    _du
end

end