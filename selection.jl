using StatsBase

function fullyrandom(individuals::Vector{Individual}, num_individuals::Integer)
    indices = sample(1:size(individuals, 1), num_individuals, replace=false)
    return individuals[indices]
end

function roulette(individuals::Vector{Individual}, num_individuals::Integer)
    weights = [individual.fitness != -Inf ? 1e-20 .+ (1 ./ (1 .+ exp.(-individual.fitness))) : 1e-20 for individual in individuals]
    nw = length(weights)
    indices = sample(1:nw, ProbabilityWeights(weights), num_individuals, replace=false)
    return individuals[indices]
end

function ranking(individuals::Vector{Individual}, num_individuals::Integer)
    weights = Vector{Float64}()
    total_individuals = length(individuals)
    selection_pressure = 1.5

    for i in 1:total_individuals
        push!(weights, 1/total_individuals * (selection_pressure - (2 * selection_pressure - 2) * ((i - 1)/(total_individuals - 1))) )
    end
    
    indices = sample(reverse(1:total_individuals), Weights(weights, sum(weights)), num_individuals, replace=false)
    return individuals[indices]
end

function tournament(individuals::Vector{Individual}, num_individuals::Int64, k)
    winners = Vector{Individual}()
    for match in 1:num_individuals
        contestants = individuals[sample(1:size(individuals, 1), k, replace=false)]
        sort!(contestants, by = v -> v.fitness, rev = true)
        push!(winners, contestants[1])
    end
    return winners
end

function byclass(individuals::Vector{Individual}, num_individuals::Int64, percent_best, percent_worst)
    if (percent_best + percent_worst) > 1.0
        throw(error("(percent_best + percent_worst) > 1.0! It must be 0 < (percent_best + percent_worst) <= 1.0"))
    end

    n_best_individuals = ceil(Int64, num_individuals * percent_best)
    n_worst_individuals = ceil(Int64, num_individuals * percent_worst)

    best_individuals = individuals[1:n_best_individuals]
    worst_individuals = individuals[(size(individuals, 1) - n_worst_individuals + 1):end]

    n_remaining_individuals = num_individuals - n_best_individuals - n_worst_individuals
    remaining_individuals = Vector{Individual}()

    println((size(individuals, 1), num_individuals, n_best_individuals, n_worst_individuals, n_remaining_individuals))
    
    if n_remaining_individuals > 0
        remaining_individuals = sample(individuals, n_remaining_individuals, replace=true) # sampling from the complete individual set with replacement.
    end

    return vcat(best_individuals, worst_individuals, remaining_individuals)
end

# It is assumed that the "individuals" is sorted by fitness, from greatest to smallest.
function select(environment::Environment, individuals::Vector{Individual}, num_individuals::Integer, method::String)
    if num_individuals > size(individuals, 1)
        throw(error("Trying to select more individuals than available."))
    end
    
    selected = Vector{Individual}()

    if method == "fullyrandom"
        selected = fullyrandom(individuals, num_individuals)
    elseif method == "roulette"
        selected = roulette(individuals, num_individuals)
    elseif method == "ranking"
        selected = ranking(individuals, num_individuals)
    elseif method == "tournament"
        selected = tournament(individuals, num_individuals, environment.selection_params[1])
    elseif method == "byclass"
        selected = byclass(individuals, num_individuals, environment.selection_params[1], environment.selection_params[2])
    else
        println("Unknow selection method ", method)
    end

    return selected
end