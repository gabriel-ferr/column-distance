##
##          Script to compute distance between a columnwise distribution.
function func_distances(dist)
    sz = size(dist)
    divergence = zeros(sz[2] - 1)
    distances = zeros(sz[2] - 1)

    for i in eachindex(divergence)
        divergence[i] = js_divergence(dist[:, i], dist[:, i + 1])
        distances[i] = sqrt(divergence[i])
    end

    points = zeros(floor(Int, length(divergence) / 50))
    divergence_means = zeros(length(points))
    distances_means = zeros(length(points))

    for i in eachindex(points)
        points[i] = 1 + (i - 1) * 50
        divergence_means[i] = mean(divergence[1 + (i - 1) * 50:i * 50])
        distances_means[i] = mean(distances[1 + (i - 1) * 50:i * 50])
    end

    return divergence, distances, points, divergence_means, distances_means
end