using Plots

function plot_solution(solution, detector_hits; vars=(:x, :z), n_particles=100, scatter_args=(), plot_args=())
    n_particles = min(n_particles, length(solution))
    plot()
    for i=1:n_particles
        plot!(solution.u[i], vars=vars, label=""; plot_args...)
    end
    for detector in detector_hits
        scatter!(
            [[hits[var] for hits in detector][begin:n_particles] for var in vars]...,
            label=""; scatter_args...
        )
    end
    plot!()
end
