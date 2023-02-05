using StatsBase

function binary_random(ga::GeneticAlgorithm)
    population = Vector{Individual}()
    
    for i in 1:ga.population_size
        chromosome = zeros(Int64, ga.n)
        indices = sample(1:ga.n, ga.s, replace=false)
        chromosome[indices] .= 1
        println(chromosome)
        individual = Individual(ga, chromosome)
        
        push!(population, individual)
    end

    return population
end

function initialize_population(ga::GeneticAlgorithm)
    population = Vector{Individual}()
    if ga.initialization_method == "binary_random"
        population = binary_random(ga)
    elseif ga.initialization_method == "something"
        population = population # placeholder
    else
        println("Unknow initialization method ", ga.initialization_method)
    end
    return population
end