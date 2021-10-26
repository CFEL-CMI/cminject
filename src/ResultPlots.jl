using Plots

function plot_solution(solution, detector_hits; vars=(:x, :z), n_particles=100)
    plot()
    for i=1:n_particles
        plot!(solution.u[i], vars=(:x, :z), label="")
    end
    for detector in detector_hits
        scatter!([[hits[var] for hits in detector][begin:n_particles] for var in vars]..., label="")
    end
    plot!()
end