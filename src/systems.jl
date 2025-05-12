##
##          Script that defines the dynamical systems used.
## -----------------------------------------------------------------------------------------------------
##      Uniform Distribution
function uniform(μ, σ, K)
    return rand(Uniform(μ, σ), K)
end
## -----------------------------------------------------------------------------------------------------
##      Normal Distribution / White Noise
function normal(μ, σ, K)
    return rand(Normal(μ, σ), K)
end
## -----------------------------------------------------------------------------------------------------
##      AR with memory rate `a`
function ar1(μ, σ, a, K)
    unif = uniform(μ, σ, K)
    norm = normal(μ, σ, K)

    result = zeros(K)
    result[1] = norm[1]

    for i in 2:K
        result[i] = (1 - a) * norm[i] + a * result[i - 1] #+ (i > 10 ? a * (unif[i - 10] / 10) * result[i - 10] : 0)
    end

    return result
end
## -----------------------------------------------------------------------------------------------------
##      Function Sin(x)
function func_sin(K::Int, ω::Float64, A::Float64)
    t = ω * π * range(0, 1.0, K)
    return A * sin.(t)
end
## -----------------------------------------------------------------------------------------------------
##      Lorenz System
function lorenz(K, σ, ρ, β; u0 = rand(3), tspan = (0.0, 20000.0), transient = 10000)
    
    function lorenz!(du, u, p, dt)
        x, y, z = u
        
        du[1] = σ * (y - x)
        du[2] = x * (ρ - z) - y
        du[3] = x * y - β * z
    end

    prob = ODEProblem(lorenz!, u0, tspan)
    sol = solve(prob)
    return prepare(sol, 0.2; transient = transient, max_length = K)
end
## -----------------------------------------------------------------------------------------------------
##      Rössler System
function rossler(K, a, b, c; u0 = rand(3), tspan = (0.0, 30000.0), transient = 10000)

    function rossler!(du, u, p, dt)
        x, y, z = u

        du[1] = - y - z
        du[2] = x + a * y
        du[3] = b + z * (x - c)
    end

    prob = ODEProblem(rossler!, u0, tspan)
    sol = solve(prob)
    return prepare(sol, 0.18; transient = transient, max_length = K)
end