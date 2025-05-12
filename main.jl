##
##          Import the necessary libraries.
using RMA
using Distances
using CairoMakie
using Statistics
using Distributions
using ProgressMeter
using RecurrenceAnalysis
using DifferentialEquations

##
##          Include some source files...
include("src/systems.jl")
include("src/distances.jl")

##
##          Configurations
data_len = 10003

##
##          Generate our data
data_lorenz = lorenz(data_len, 10.0, 28.0, 8.0/3.0)
data_rossler = rossler(data_len, 0.2,  0.2, 18.0)

data = [
    [
        "Uniform Distribution",
        [1, 1],
        uniform(0, 1, data_len)
    ],
    [
        "Normal Distribution",
        [1, 2],
        normal(0, 1, data_len)
    ],
    [
        "Random with Memory",
        [1, 3],
        ar1(0, 1, 0.8, data_len)
    ],
    [
        "Function sin(4.2π t)",
        [2, 1],
        func_sin(data_len, 4.2, 1.0)
    ],
    [
        "Lorenz System",
        [2, 2],
        data_lorenz
    ],
    [
        "Rössler System",
        [2, 3],
        data_rossler
    ],
    [
        "Lorenz System component x",
        [3, 1],
        data_lorenz[1, :]
    ],
    [
        "Lorenz System component y",
        [3, 2],
        data_lorenz[2, :]
    ],
    [
        "Lorenz System component z",
        [3, 3],
        data_lorenz[3, :]
    ],
    [
        "Rössler System component x",
        [4, 1],
        data_rossler[1, :]
    ],
    [
        "Rössler System component y",
        [4, 2],
        data_rossler[2, :]
    ],
    [
        "Rössler System component z",
        [4, 3],
        data_rossler[3, :]
    ]
]

data_lorenz = Nothing
data_rossler = Nothing

##      Compute our graphics
@showprogress for n = 2:4
    fig_div = Figure(size = (1800, 2400))
    fig_met = Figure(size = (1800, 2400))

    counting = 0
    @showprogress for d in data
        title, pos, set = d

        ##      Get the th that maximize recurrence entropy.
        th, _ = find_parameters(set, 2; threshold_max = maximum(pairwise(Euclidean(), set, set)))
        dist = distribution(set, th, n; sampling_mode = :columnwise_full)

        divergence, distances, points, divergence_means, distances_means = func_distances(dist)
        dist = Nothing

        ax_div = Axis(fig_div[pos...], title = title, titlealign = :left, limits = (nothing, (0.00001, 1.0)), xlabel = "Time", ylabel = "log[JS Divergence]", yscale = log)
        lines!(ax_div, divergence, color = :tomato, linewidth = 1)
        scatterlines!(ax_div, points, divergence_means, color = :black, markersize = 4, linewidth = 2)

        ax_met = Axis(fig_met[pos...], title = title, titlealign = :left, limits = (nothing, (0.00001, 1.0)), xlabel = "Time", ylabel = "log[sqrt[JS Divergence]] ", yscale = log)
        lines!(ax_met, distances, color = :tomato, linewidth = 1)
        scatterlines!(ax_met, points, distances_means, color = :black, markersize = 4, linewidth = 2)

        ##      I will do some other grapics too...
        if (pos == [2, 2] || pos == [2, 3])
            println(string("\nIgnoring: ", title))
            continue
        end

        fig_serie = Figure(size = (600, 600))

        ax_res = Axis(fig_serie[1, 1], title = title, titlealign = :left, limits = (nothing, (0.00001, 1.0)), ylabel = "log[sqrt[JS Divergence]] ", yscale = log)
        hidexdecorations!(ax_res)
        lines!(ax_res, distances, color = :tomato, linewidth = 1)
        scatterlines!(ax_res, points, distances_means, color = :black, markersize = 4, linewidth = 2)

        ax_serie = Axis(fig_serie[2, 1], title = title, titlealign = :left, ylabel = "value", xlabel = "Time")
        lines!(ax_serie, set, color = :blue, linewidth = 1)

        rowsize!(fig_serie.layout, 2, Relative(0.3))

        save(string("img/n_", n, "_serie_", counting, ".png"), fig_serie)
        counting += 1
    end

    save(string("img/n_", n, "_div.png"), fig_div)
    save(string("img/n_", n, "_met.png"), fig_met)

    println("\nDone !")
end