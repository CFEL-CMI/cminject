using CMInject
using ArgParse
using Distributions
using Plots
using DifferentialEquations

#=
plotly()
sol, detector_hits = CMInject.ensemble_solve(100)
#plot = CMInject.plot_solution(sol, detector_hits; vars=(:x, :z, :vz))
display(plot)
=#


# Wrap the Distribution in a new struct, so we won't accidentally come in conflict with some other
# implementation of ArgParse.parse_item(::Type{Distribution}, x::AbstractString)
struct ArgParseDistribution{Dist<:Distribution}
    dist::Dist
end

function parse_distribution(x::AbstractString) where Dist<:Distribution
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
        else
            error("I don't know this distribution type: $(x[1])")
            # TODO implement more distributions with some elseif cases...
        end
    end
end

# Shorthand for extracting the contained Distribution from an ArgParseDistribution
_d(x::ArgParseDistribution) = x.dist

function ArgParse.parse_item(::Type{ADist}, x::AbstractString) where ADist<:ArgParseDistribution
    return ArgParseDistribution(parse_distribution(x))
end

function run_argparse(args)
    s = ArgParseSettings()
    @add_arg_table! s begin
        "-f"
            arg_type = String
            help = "Path of a flow-field HDF5 file"
            required = true
        "-n"
            arg_type = Int
            help = "The number of particles to simulate"
            required = true
        "--x"
            arg_type = ArgParseDistribution
            help = "x (position) distribution"
            required = true
        "--z"
            arg_type = ArgParseDistribution
            help = "z (position) distribution"
            required = true
        "--vx"
            arg_type = ArgParseDistribution
            help = "vx (velocity) distribution"
            required = true
        "--vz"
            arg_type = ArgParseDistribution
            help = "vz (velocity) distribution"
            required = true
        "--r"
            arg_type = ArgParseDistribution
            help = "Radius distribution"
            required = true
        "--rho"
            arg_type = ArgParseDistribution
            help = "Material density distribution"
            required = true
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
            help = "Add this option to plot 100 of the simulated particles and detector hits"
            action = :store_true
    end
    return parse_args(args, s)
end

function main()
    args = run_argparse(ARGS)

    field = StokesFlowField(args["f"], args["fT"], args["fM"], args["fMu"])

    dists = (
        x  = _d(args["x"]),  z = _d(args["z"]),
        vx = _d(args["x"]), vz = _d(args["vz"]),
        r  = _d(args["r"]),  ρ = _d(args["rho"])
    )
    ParticleType = SphericalParticle2D{Float64}
    source = SamplingSource{ParticleType}(dists)

    experiment = Experiment(;
        source=source,
        n_particles=args["n"],
        fields=(field,),
        detectors=Tuple(SectionDetector{Float64,:z}.(args["d"], true)),
        time_span=Tuple(args["t"]), time_step=args["dt"],
        solver=EulerHeun(), ensemble_alg=EnsembleThreads()
    )
    solution, hits, particles = simulate(experiment)

    if args["plot"]
        plotly()
        display(plot_results(solution, hits))
    end

    if !isempty(args["o"])
        store = HDF5ResultStorage{ParticleType}(args["o"])
        store_results!(store, solution, hits, particles; store_trajectories=args["T"])
    end
end

main()