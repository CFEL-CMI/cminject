using DifferentialEquations
import SciMLBase: AbstractSDEAlgorithm, EnsembleAlgorithm
using Parameters

"""
    Experiment(; source, n_particles, time_span, time_step
               [fields=(), detectors=(), solver=EulerHeun(),
                solver_opts=(adaptive=false, dense=false, unstable_check=always_false),
                ensemble_alg=EnsembleThreads()])

Fully specifies a particle trajectory simulation experiment in CMInject, consisting of:

- `source`, a source of particles
- `n_particles`, the number of particles to be simulated
- `time_span`, the timespan to simulate over
- `time_step`, the macro-timestep for the integrator
- `fields`, a tuple of (acceleration) fields - none by default
- `detectors`, a tuple of detectors - none by default
- `solver`: A solver algorithm for `SDEProblem`s
- `solver_opts`: A (Named)Tuple of options to be passed to the solver
- `ensemble_alg`: The ensemble algorithm to use for parallelizing over the particles,
   see the DifferentialEquations.jl docs
"""
@with_kw struct Experiment{
    Src<:AbstractSource,Fields,Dets,
    TimeStep,SolvAlg<:AbstractSDEAlgorithm,SolvOpts,EnsAlg<:EnsembleAlgorithm
}
    source::Src
    n_particles::Int64
    time_span::Tuple{TimeStep,TimeStep}
    time_step::TimeStep
    fields::Fields = ()
    detectors::Dets = ()

    solver::SolvAlg = EulerHeun()
    solver_opts::SolvOpts = (adaptive=false, dense=false, unstable_check=always_false)
    ensemble_alg::EnsAlg = EnsembleThreads()
end


"""
    f!(du, u, params, t)

The deterministic part of an SDEProblem, as a generated function which
specializes on a particular particle type and collection of fields.
"""
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
"""
    g!(du, u, params, t)

The stochastic part of an SDEProblem, as a generated function which
specializes on a particular particle type and collection of fields.

!!! note "Support of noise types"
    Currently only supports diagonal additive noise, and will sum
    the noise terms of all locally applicable fields together.
"""
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

"""
    always_false(args...)

Ignores its arguments and always returns false. Useful when a function returning a boolean value is expected and
you always want to return false. CMInject currently uses it to pass as an `unstable_check` as we are currently
doing NaN-based handling.
"""
@inline function always_false(args...)
    false
end

"""
    make_initial_values(source, fields, n::Integer)

From a `source` and a collection of `fields`, generates `n` initial particle values and returns a tuple of
`(u₀, p₀, prob_func)`, where

- `u₀` is the phase-space position of the first generated particle
- `p₀` is the set of particle properties of the first generated particle
- `prob_func` is a function that returns an updated DEProblem within an ensemble, as required by `EnsembleProblem`.
"""
function make_initial_values(source, fields, n::Integer)
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

    u0, params, particles, prob_func
end

"""
    simulate(exp::Experiment)

Simulates an `Experiment`, returning a tuple of `(trajectories, detector_hits)`.
"""
function simulate(exp::Experiment)
    u0, params, particles, prob_func = make_initial_values(exp.source, exp.fields, exp.n_particles)
    prob = SDEProblem{true}(f!, g!, u0, exp.time_span, params)
    ens_prob = EnsembleProblem(prob, prob_func=prob_func, safetycopy=false)

    solution = solve(
        ens_prob, exp.solver, exp.ensemble_alg;
        trajectories=exp.n_particles,
        dt=exp.time_step,
        exp.solver_opts...
    )
    hits = [calculate_hits(solution, detector) for detector in exp.detectors]
    return solution, hits, particles
end
