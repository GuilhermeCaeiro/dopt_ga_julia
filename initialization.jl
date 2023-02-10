using StatsBase

function binary_random(environment::Environment)
    population = Vector{Individual}()
    
    for i in 1:environment.population_size
        chromosome = zeros(Int64, environment.n)
        indices = sample(1:environment.n, environment.s, replace=false)
        chromosome[indices] .= 1
        #println(chromosome)
        individual = Individual(environment, chromosome)
        
        push!(population, individual)
    end

    return population
end

function binary_biased(environment::Environment, min_diff::Float64)
    chromosomes = Vector{Vector{Int64}}()
    
    while size(chromosomes, 1) < environment.population_size
        chromosome = zeros(Int64, environment.n)
        chromosome[sample(1:size(chromosome. 1), environment.s, replace = false)] .= 1
        found_close = false

        for existing_chromosome in chromosomes
            if mean(abs(existing_chromosome - chromosome)) < min_diff
                found_close = true
                break
            end
        end
        if !found_close
            push!(chromosomes, chromosome)
        end
    end

    return chromosomes
end




function initialize_population(environment::Environment)
    population = Vector{Individual}()
    if environment.initialization_method == "binary_random"
        population = binary_random(environment)
    elseif environment.initialization_method == "something"
        population = population # placeholder
    else
        println("Unknow initialization method ", environment.initialization_method)
    end
    return population
end