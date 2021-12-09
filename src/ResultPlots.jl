using Plots

function plot_results(
    solution, detector_hits; vars=(:x, :z), n_particles=100,
    plot_detectors=true, plot_trajectories=true, scatter_args=(), plot_args=()
)
    plot()

    if plot_trajectories
        n = min(n_particles, length(solution))
        for i=1:n
            plot!(solution.u[i], vars=vars, label=""; plot_args...)
        end
    end
    if plot_detectors
        for detector in detector_hits
            n = min(n_particles, length(detector))
            scatter!(
                [[hits[var] for hits in detector][begin:n] for var in vars]...,
                label=""; scatter_args...
            )
        end
    end

    plot!()
end
