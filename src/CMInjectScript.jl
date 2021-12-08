## NOTE: This is not a functional part of CMInject.jl. It's an initial development notebook, to learn about
## performance issues and implement a working ALS simulation by hand, with no generic framework in sight.
## I'm keeping it around as long as there are features to be carried over to the framework. -Simon
##
## Features to transfer:
##  - 3D simulations (in particular: with an underlying axisymmetric 2D field)
##  - Detectors (post-processing style like done here)
##  - Some plotting utilities

module CMInjectScript

using HDF5
using DifferentialEquations
using StaticArrays
using DiffEqGPU
using Plots
using Statistics

const kB = Float64(1.380649e-23)

# D = spatial dimensions, N = interpolated quantities
struct RegularGrid{D,N,T,A <: AbstractArray{T, D}}
    ranges::NTuple{D, StepRangeLen{T, T, T}}
    data::NTuple{N, A}
end

function interpol(grid::RegularGrid{2,N,T}, xs::V, ys::V) where {N, T, M<:AbstractMatrix{T}, V<:AbstractVector{T}}
    out = Array{T}(undef, length(xs), N)
    interpol!(out, grid, xs, ys)
    out
end

function interpol!(out::M, grid::RegularGrid{2,N,T}, xs::V, ys::V) where {N, T, M<:AbstractMatrix{T}, V<:AbstractVector{T}}
    for (k, (x, y)) in enumerate(zip(xs, ys))
        outv = @view out[k,:]
        interpol1!(outv, grid, x, y)
    end
    nothing
end

function interpol1!(out::V, grid::RegularGrid{2,N,T}, x::T, y::T) where {N, T, V<:AbstractVector{T}}
    xrange = grid.ranges[1]
    yrange = grid.ranges[2]

    if isnan(x) || isnan(y)
        for i in 1:N
            out[i] = NaN
        end

        return nothing
    end

    scaledx = (x - xrange[1]) / xrange.step
    scaledy = (y - yrange[1]) / yrange.step
    xw, xi = modf(scaledx)
    yw, yi = modf(scaledy)
    xi = Int(xi)+1  # TODO: +1 correct? seems good at first glance, also in 2D, but hmm...
    yi = Int(yi)+1

    outofbounds = (xi < 1) || (xi >= xrange.len) || (yi < 1) || (yi >= yrange.len)

    if outofbounds
        for i in 1:N
            out[i] = NaN
        end
    else
        for i in 1:N
            dat = grid.data[i]'
            a0 = dat[xi,yi  ]*(1-xw) + dat[xi+1,yi  ]*xw
            a1 = dat[xi,yi+1]*(1-xw) + dat[xi+1,yi+1]*xw
            out[i] = a0*(1-yw) + a1*yw
        end
    end

    nothing
end

function to_reused_interpol1(grid::RegularGrid{D,N,T}) where {D,N,T}
    vec = collect(repeat([NaN::T], N))  # TODO use MVector -- but the GPU code doesn't like that
    return (x::T, y::T) -> begin
        interpol1!(vec, grid, x, y)
        vec
    end
end

function to_reused_interpol(grid::RegularGrid{D,N,T}, M::Integer) where {D,N,T}
    matrix = Matrix{T}(undef, M, N)
    return (xs::Vector{T}, ys::Vector{T}) -> begin
	interpol!(matrix, grid, xs, ys)
	matrix
    end
end

function nodes_to_range(nodes::AbstractArray{T,N}; delta=nothing) where {T, N}
    xα = nodes[1]
    xβ = nodes[2]
	Δx = isnothing(delta) ? xβ - xα : delta
    StepRangeLen{T, T, T}(xα, Δx, length(nodes))
end

function hdf5_to_regulargrid(filename::AbstractString)
    h5open(filename) do h5f
	idx = h5f["index"]

	index = [read(idx, k) for k in keys(idx)]
	index = tuple(index...)
	sizes = (size(idx)[1] for idx in reverse(index))

	dat = h5f["data"]
	data = [
	    permutedims(reshape(read(dat, k), sizes...), (2,1))'
	    for k in keys(dat)
	]

	ranges = nodes_to_range.(index)
	RegularGrid(ranges, Tuple(data))
    end
end

function stokes!(du, u, params, t)
    x, vx, z, vz = u
    itp, rₚ, ρ, μ, T, mᶠ = params

    fluid = itp(x, z)
    vxᶠ, vzᶠ, pᶠ = fluid
    Δvx = vxᶠ - vx
    Δvz = vzᶠ - vz

    Knf = μ * √(π * kB * T / (2mᶠ))
    mp = 4/3 * rₚ^3. * π * ρ
    Kn = Knf / (pᶠ * rₚ)
    Cc = 1 + Kn * (1.231 + 0.4695 * exp(-1.1783 / Kn))
    a₀ = 6π * μ * rₚ / Cc / mp

    @inbounds begin
        du[1] = u[2]
        du[2] = a₀ * Δvx
        du[3] = u[4]
        du[4] = a₀ * Δvz
    end
    nothing
end

function ahmadi!(du, u, params, t)
    x, vx, z, vz = u
    itp, rₚ, ρ, μ, T, mᶠ = params

    fluid = itp(x, z)
    _, _, pᶠ = fluid
    Knf = μ * √(π * kB * T / (2mᶠ))
    Kn = Knf / (pᶠ * rₚ)
    Cc = 1 + Kn * (1.231 + 0.4695 * exp(-1.1783 / Kn))
    s₀fπ = π * 216μ * kB * T / (π^2 * (2rₚ)^5 * ρ^2)
    s₀ = √(s₀fπ / Cc)

    @inbounds begin
        du[1] = 0
        du[2] = s₀
        du[3] = 0
        du[4] = s₀
    end
    nothing
end

function stokes3!(du, u, params, t)
    x, vx, y, vy, z, vz = u
    itp, rₚ, ρ, μ, T, mᶠ = params

    fluid = itp(√(x^2 + y^2), z)
    vrᶠ, vzᶠ, pᶠ = fluid

    if isnan(pᶠ)
        @inbounds begin
            du[1] = u[2]
            du[2] = 0
            du[3] = u[4]
            du[4] = 0
            du[5] = u[6]
            du[6] = 0
        end
    else
        ∠xy = atan(y, x)
        vxᶠ = cos(∠xy) * vrᶠ
        vyᶠ = sin(∠xy) * vrᶠ
        Δvx = vxᶠ - vx
        Δvy = vyᶠ - vy
        Δvz = vzᶠ - vz

        Knf = μ * √(π * kB * T / (2mᶠ))
        mp = 4/3 * rₚ^3 * π * ρ
        Kn = Knf / (pᶠ * rₚ)
        Cc = 1 + Kn * (1.231 + 0.4695 * exp(-1.1783 / Kn))
        a₀ = 6π * μ * rₚ / Cc / mp
        @inbounds begin
            du[1] = u[2]
            du[2] = a₀ * Δvx
            du[3] = u[4]
            du[4] = a₀ * Δvy
            du[5] = u[6]
            du[6] = a₀ * Δvz
        end
    end

    nothing
end

function ahmadi3!(du, u, params, t)
    x, vx, y, vy, z, vz = u
    itp, rₚ, ρ, μ, T, mᶠ = params

    fluid = itp(√(x^2 + y^2), z)
    _, _, pᶠ = fluid
    Knf = μ * √(π * kB * T / (2mᶠ))
    Kn = Knf / (pᶠ * rₚ)
    Cc = 1 + Kn * (1.231 + 0.4695 * exp(-1.1783 / Kn))
    s₀fπ = π * 216μ * kB * T / (π^2 * (2rₚ)^5 * ρ^2)
    s₀ = √(s₀fπ / Cc)

    if isnan(s₀)
        s₀ = 0.0
    end

    @inbounds begin
        du[1] = 0
        du[2] = s₀
        du[3] = 0
        du[4] = s₀
        du[5] = 0
        du[6] = s₀
    end
    nothing
end

const example_grid = hdf5_to_regulargrid("/home/chip/Documents/cminject-materials/goldspheres-27nm/nitrogen_1.8mbar_extended_nan.h5")
const example_ipl = to_reused_interpol1(example_grid)

const μ = 1.76e-5
const rp = 13.5e-9
const ρ = 19320.0
const T = 293.15
const mf = 4.27e-26

const example_params = (example_ipl, rp, ρ, μ, T, mf)

const example_u0 = [0.0, 0.0, -0.128, 10.0]
const example_problem = SDEProblem{true}(stokes!, ahmadi!, example_u0, (0.0, 0.03), example_params)
const example_ode_problem = ODEProblem{true}(stokes!, example_u0, (0.0, 0.03), example_params)
const example_ensemble_problem = EnsembleProblem(
    example_problem,
    prob_func=(prob,i,repeat) -> begin
        u0 = [randn()*0.001*1.41, randn()*1.41, -0.128, 10.0+randn()]
        params = (to_reused_interpol1(example_grid), rp, ρ, μ, T, mf)
        SDEProblem{true}(stokes!, ahmadi!, u0, (0.0, 0.03), params)
    end,
    safetycopy=true
)

const example_u0_3d = [0.0, 0.0, 0.0, 0.0, -0.128, 10.0]
const example_problem_3d = SDEProblem{true}(stokes3!, ahmadi3!, example_u0_3d, (0.0, 0.1), example_params)
const example_ensemble_problem_3d = EnsembleProblem(
    example_problem_3d,
    prob_func=(prob,i,repeat) -> begin
        u0 = [randn()*0.001, randn(), randn()*0.001, randn(), -0.128, 10.0+randn()]
        params = (to_reused_interpol1(example_grid), rp, ρ, μ, T, mf)
        SDEProblem{true}(stokes3!, ahmadi3!, u0, (0.0, 0.03), params)
    end,
    safetycopy=false
)

const always_false = (_...) -> false

function example_solve()
    solve(example_problem, EulerHeun(), dt=1e-5, adaptive=false, dense=false, unstable_check=always_false)
end

function example_solve_ensemble(n; ensemblealg=EnsembleThreads(), adaptive=false, dense=false)
    solve(example_ensemble_problem, EulerHeun(), ensemblealg,
          dt=1e-5, adaptive=adaptive, dense=dense, unstable_check=always_false;
          trajectories=n)
end

function example_solve_ode()
    solve(example_ode_problem, Heun(), dt=1e-5, adaptive=false, dense=false, unstable_check=always_false)
end

iscallable(f) = !isempty(methods(f))

function plot_solutions!(solutions, x, y; z=nothing, k=length(solutions), lastn=nothing)
    #plot(size=size)
    for i=1:k
        U = hcat(solutions[i].u...)
        nanidx = findfirst(isnan, U[1,:])
        if nanidx !== nothing
            U = U[:, begin:nanidx]
        end

        if lastn !== nothing
            U = U[end-lastn:end, :]
        end

        x_ = iscallable(x) ? x(U) : U[x,:]
        y_ = iscallable(y) ? y(U) : U[y,:]
        if z !== nothing
            z_ = iscallable(z) ? z(U) : U[z,:]
            plot!(x_, y_, z_)
        else
            plot!(x_, y_)
        end
    end
end

function get_detectors(solutions, zs; mode=:d3)
    zidx = mode == :d3 ? 5 : 3

    detectors = []
    for z in zs
        detector = hcat(get_at_zpos.(solutions[:], z; zidx=zidx)...)
        push!(detectors, detector)
    end
    detectors
end


function plot_detectors!(detectors; mode=:d3)
    for det in detectors
        if mode == :d3
            show = (det[1, :], det[3, :], det[5, :])
        else
            show = (det[1, :], det[3, :])
        end
        scatter!(show...; markersize=1)
    end
end


function get_at_zpos(trajectory, zpos; zidx=5)
    z = trajectory[zidx, :]
    i₁ = findfirst(z .>= zpos)
    if isnothing(i₁)
        return zeros(Float64, size(trajectory)) / 0
    end

    # there are some unhandled edge cases here...
    i₀ = i₁ - 1

    w₁ = (zpos - z[i₀]) / (z[i₁] - z[i₀])
    w₀ = 1 - w₁

    w₀ * trajectory[:, i₀] .+ w₁ * trajectory[:, i₁]
end

function get_detector(trajectories, zpos; zidx=5)
    hcat(get_at_zpos.(trajectories[:], zpos; zidx=zidx)...)
end

function get_detector_r(detector)
    sqrt.(detector[1, :].^2 .+ detector[3, :].^2)
end

function get_detector_x(detector)
    detector[1, :]
end


function main(n=1000)
    plotly()

    myplot = plot()
    sol = solve(example_ensemble_problem_3d, EulerHeun(), dt=1e-5,
                dense=false, unstable_check=always_false; trajectories=n)
    plot_solutions!(sol, 1, 3, z=5, k=250)

    zs = -5e-3:1e-3:5e-3
    detectors = get_detectors(sol, zs, mode=:d3)
    plot_detectors!(detectors, mode=:d3)

    #myplot

    rs = [filter(!isnan, get_detector_r(det)) for det in detectors]
    focus_sizes = 4*Statistics.quantile.(rs, 0.7)
    display(plot(zs, focus_sizes))
end

end