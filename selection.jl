using StatsBase

function fullyrandom(individuals::Vector{Individual}, num_individuals::Integer)
    #println(size(individuals), " ", num_individuals)
    indices = sample(1:size(individuals, 1), num_individuals, replace=false)
    return individuals[indices]
end

#
function roulette(individuals::Vector{Individual}, num_individuals::Integer)
    weights = Vector{Integer}()

    for i in 1:size(individuals, 1)
        push!(weights, individuals[i].fitness)
    end

    weights = weights / sum(weights)
    
    indices = sample(1:size(individuals, 1), Weights(weights), num_individuals, replace=false)
    return individuals[indices]
end

# TODO improve
function ranking(individuals::Vector{Individual}, num_individuals::Integer)
    weights = Vector{Integer}()

    for i in reverse(1:size(individuals, 1))
        push!(weights, i)
    end
    
    indices = sample(1:size(individuals, 1), Weights(weights), num_individuals, replace=false)
    return individuals[indices]
end

# It is assumed that the "individuals" is sorted by fitness, from greatest to smallest.
function select(ga::GeneticAlgorithm, individuals::Vector{Individual}, num_individuals::Integer, method::String)
    selected = Vector{Individual}()

    if method == "fullyrandom"
        selected = fullyrandom(individuals, num_individuals)
    elseif method == "roulette"
        selected = roulette(individuals, num_individuals)
    elseif method == "ranking"
        selected = ranking(individuals, num_individuals)
    else
        println("Unknow selection method ", method)
    end

    return selected
end