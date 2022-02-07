using CMInject
using ArgParse
using Distributions
using Plots
using DifferentialEquations


# Wrap the Distribution in a new struct, so we won't accidentally come in conflict with some other
# implementation of ArgParse.parse_item(::Type{Distribution}, x::AbstractString)
struct ArgParseDistribution{Dist<:Distribution}
    dist::Dist
end

"""
    parse_distribution(x::AbstractString)

Parses a `Distribution` from an `AbstractString` and returns it.
Calls `error` if parsing is unsuccessful. Currently supports the following:

    - "G[0.0,1.0]"  -- a Gaussian/Normal distribution with μ=0.0, σ=1.0
    - "D[0,2]"      -- a Discrete Uniform distribution from between a=0 and b=2
    - "-0.1285"     -- a Dirac distribution always sampling to `-0.1285`

Spaces within the expressions are not permitted.

!!! note
    All returned continuous (and Dirac) distributions use Float64 values.
"""
function parse_distribution(x::AbstractString)
    if !isletter(x[1])
        # Probably a number, try to parse it as one, then return a Dirac distribution
        num = parse(Float64, x)
        return Dirac(num)
    else
        # We seem to have something like G[0,1] on our hands...
        # check that the arguments are neatly wrapped in [] brackets
        if !(x[2] == '[' && x[end] == ']')
            error("I expected a distribution specification like G[0,1] but you passed $(x)")
        end

        # Pull out the bit of the string that should be the distribution arguments
        args = x[3:end-1]
        if x[1] == 'G'
            # Gaussian/Normal distribution
            mu, sigma = parse.(Float64, split(args, ','))
            return Normal(mu, sigma)
        elseif x[1] == 'D'
            # Discrete uniform distribution
            a, b = parse.(Int32, split(args, ','))
            return DiscreteUniform(a, b)
        else
            error("I don't know this distribution type: $(x[1])")
            # TODO implement more distributions with some elseif cases...
        end
    end
end

"""
    _d(x::ArgParseDistribution{Dist})::Dist = x.dist

Shorthand for extracting the contained Distribution from an ArgParseDistribution
"""
function _d(x::ArgParseDistribution{Dist})::Dist where {Dist<:Distribution}
    x.dist
end

function ArgParse.parse_item(::Type{ADist}, x::AbstractString) where ADist<:ArgParseDistribution
    return ArgParseDistribution(parse_distribution(x))
end

"""
    run_argparse(args)

Defines an ArgParseSettings for this script and runs it on the input argument args.
For parsing from the command-line, pass the global variable ARGS.
"""
function run_argparse(args)
    s = ArgParseSettings()
    @add_arg_table! s begin
        "-f"
            arg_type = String
            help = "Path of a flow-field HDF5 file"
            default = ""
        "-e"
            arg_type = String
            help = "Path of an electric-field HDF5 file"
            default = ""
        "-s"
            arg_type = String
            help = "Path of a Stark-curve HDF5 file"
            default = ""
        "-n"
            arg_type = Int
            help = "The number of particles to simulate"
            required = true
        "--dn"
            arg_type = Int
            help = "The number of dimensions"
            default = 2
        "--x"
            arg_type = ArgParseDistribution
            help = "x (position) distribution"
            required = true
        "--y"
            arg_type = ArgParseDistribution
            help = "y (position) distribution. Only necessary for 3 dimensions"
            default = ArgParseDistribution(Dirac(0))
        "--z"
            arg_type = ArgParseDistribution
            help = "z (position) distribution"
            required = true
        "--vx"
            arg_type = ArgParseDistribution
            help = "vx (velocity) distribution"
            required = true
        "--vy"
            arg_type = ArgParseDistribution
            help = "vy (velocity) distribution. Only necessary for 3 dimensions"
            default = ArgParseDistribution(Dirac(0))
        "--vz"
            arg_type = ArgParseDistribution
            help = "vz (velocity) distribution"
            required = true
        "--r"
            arg_type = ArgParseDistribution
            help = "Radius distribution"
            default = ArgParseDistribution(Dirac(0))
        "--rho"
            arg_type = ArgParseDistribution
            help = "Material density distribution"
            default = ArgParseDistribution(Dirac(0))
        "--m"
            arg_type = ArgParseDistribution
            help = "Mass of the particle. Use this parameter for Stark effect"
            default = ArgParseDistribution(Dirac(0))
        "--J"
            arg_type = ArgParseDistribution
            help = "The distribution of the J quantum number"
            default = ArgParseDistribution(CMInject.DiscreteUniform(0,0))
        "--M"
            arg_type = ArgParseDistribution
            help = "The distribution of the M quantum number"
            default = ArgParseDistribution(CMInject.DiscreteUniform(0,0))
        "-t"
            arg_type = Float64
            help = "The time-span (beginning and end) to simulate for"
            nargs = 2
            required = true
        "--dt"
            arg_type = Float64
            help = "The integrator time-step to use"
            default = 1e-5
        "-d"
            arg_type = Float64
            nargs = '*'
            help = "Detector z positions"
        "-o"
            arg_type = String
            help = "Output HDF5 file path"
        "-T"
            help = "Add this option to store trajectories in the output file (-o)"
            action = :store_true
        "--fT"
            arg_type = Float64
            help = "Temperature of the flow-field"
            default = 293.15
        "--fM"
            arg_type = Float64
            help = "Mass of a single gas/fluid particle of the flow-field"
            default = 4.27e-26  # N2 (nitrogen)
        "--fMu"
            arg_type = Float64
            help = "The dynamic viscosity of the flow-field"
            default = 1.76e-5  # N2 (nitrogen) at/around room temp
        "--plot"
            help = "Add this option to plot 100 simulated trajectories and detector hits"
            action = :store_true
    end
    return parse_args(args, s)
end

"""
    main()

Declares and simulates an `Experiment` with a single `StokesFlowField`
and `SphericalParticle2D{Float64}` particles, parameterized by the
command-line arguments (`ARGS`).

Returns a tuple `(solution, detectors, particles, theplot)`.
`theplot` may be `nothing` if `--plot` was not passed via ARGS.
"""
function main()
    args = run_argparse(ARGS)

    dimensions = args["dn"]
    field = Nothing
    dists = (
             x  = _d(args["x"]),  z = _d(args["z"]),
             vx = _d(args["x"]), vz = _d(args["vz"])
            )
    stateDists = (
                  J = _d(args["J"]), M = _d(args["M"])
                 )
    if (dimensions == 3)
        dists = merge(dists, (y = _d(args["y"]), vy = _d(args["vy"])))
    end

    # TODO: Merge those sections by generalizing stuff
    if (length(args["f"]) != 0 && length(args["e"]) == 0)
        # FLOW field
        # Fow now we use r and ρ for the flow field
        
        dists = merge(dists, (r = _d(args["r"]), ρ = _d(args["rho"])))
        field = StokesFlowField(args["f"], args["fT"], args["fM"], args["fMu"])
        ParticleType = dimensions == 2 ? SphericalParticle2D{Float64} : SphericalParticle3D{Float64}
        source = SamplingSource{ParticleType}(dists)

    elseif (length(args["f"]) == 0 && length(args["e"]) != 0)
        # STARK effect
        # We use m for the stark effect

        dists = merge(dists, (m = _d(args["m"]),))
        field = ElectricField(args["e"])
        ParticleType = dimensions == 2 ? StarkParticle2D{Float64} : StarkParticle{Float64}

        dictDists = Dict(pairs(dists))
        dictStateDists = Dict(pairs(stateDists))

        source = StarkSamplingSource{ParticleType, Float64}(dictDists, dictStateDists, args["s"])

    else
        error("Exactly one of the arguments of f or e has to be specified -",
              " no field or multiple fields aren't supported (yet)")
    end

    experiment = Experiment(;
        source=source,
        n_particles=args["n"],
        fields=(field,),
        detectors=Tuple(SectionDetector{Float64,:z}.(args["d"], true)),
        time_span=Tuple(args["t"]), time_step=args["dt"],
        solver=EulerHeun(), ensemble_alg=EnsembleThreads()
    )
    solution, detectors, particles = simulate(experiment)
    theplot = nothing

    if args["plot"]
        gr()
        display(plot_results(solution, detectors))

        # required so the plot doesn't immediately close... there may be better ways.
        print("Press ENTER to exit.")
        readline()
    end

    if !isnothing(args["o"])
        store = HDF5ResultStorage{ParticleType}(args["o"])
        store_results!(store, solution, detectors, particles; store_trajectories=args["T"])
    end

    return solution, detectors, particles, theplot
end

# this puts the four result variables into the global namespace when running interactively
solution, detectors, particles, theplot = main()
