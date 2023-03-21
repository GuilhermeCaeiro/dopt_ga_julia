using StatsBase

function binary_random(environment::Environment)
    population = Vector{Individual}()
    
    for i in 1:environment.population_size
        chromosome = zeros(Int64, environment.n)
        indices = sample(1:environment.n, environment.s, replace=false)
        chromosome[indices] .= 1
        individual = Individual(environment, chromosome)
        
        push!(population, individual)
    end

    return population
end

function binary_biased(environment::Environment, min_diff)
    chromosomes = Vector{Vector{Int64}}()
    
    while size(chromosomes, 1) < environment.population_size
        chromosome = zeros(Int64, environment.n)
        chromosome[sample(1:size(chromosome, 1), environment.s, replace = false)] .= 1
        found_close = false

        for existing_chromosome in chromosomes
            if mean(abs.(existing_chromosome - chromosome)) < min_diff
                found_close = true
                break
            end
        end
        if !found_close
            push!(chromosomes, chromosome)
        end
    end

    return [Individual(environment, chromosome) for chromosome in chromosomes]
end

function binary_biasedweighted(environment::Environment)
    println("binary_biasedweighted")
    chromosomes = Vector{Vector{Int64}}()
    
    chromosome = zeros(Int64, environment.n)
    chromosome[sample(1:size(chromosome, 1), environment.s, replace = false)] .= 1
    push!(chromosomes, chromosome)

    while size(chromosomes, 1) < environment.population_size
        weights = 1 ./ (sum(chromosomes) .+ 1) # the sum() is by default row-wise, producing a Vector.
        weights = weights ./ sum(weights)

        chromosome = zeros(Int64, environment.n)
        chromosome[sample(1:size(chromosome, 1), Weights(weights), environment.s, replace = false)] .= 1
        push!(chromosomes, chromosome)
    end

    return [Individual(environment, chromosome) for chromosome in chromosomes]
end

function binary_guided(environment::Environment)

    R = map(Int64, environment.R)
    m, n = size(environment.A)
    s = environment.s
    U, = svd(environment.A; full=true)
    x_bar = vec(sum(U .^ 2; dims = 2))   
    x_bar[R] .= 0.0
    chromosomes = Vector{Vector{Int64}}(undef, environment.population_size)
    x = zeros(m)
    x[R] .= 1.0 # choose independent rows
    phi = partialsortperm(x_bar, (1:s-n); rev = true)
    x[phi] .= 1.0
    zeros_idx = findall(x .== 0.0)
    ones_idx = findall(x .== 1.0)
    changes = Int(floor(m / 10))
    chromosomes[1] = deepcopy(x)
    for i in 2:environment.population_size
        _x = deepcopy(x)
        _x[sample(zeros_idx, changes, replace = false)] .= 1.0
        _x[sample(ones_idx, changes, replace = false)] .= 0.0
        chromosomes[i] = deepcopy(_x)
    end

    return [Individual(environment, chromosome) for chromosome in chromosomes]
end


function initialize_population(environment::Environment)
    population = Vector{Individual}()
    if environment.initialization_method == "binary_random"
        population = binary_random(environment)
    elseif environment.initialization_method == "binary_biased"
        # initialization_params is assumed to be  a float
        population = binary_biased(environment, environment.initialization_params[1])
    elseif environment.initialization_method == "binary_biasedweighted"
        # initialization_params is assumed to be  a float
        population = binary_biasedweighted(environment)
    elseif environment.initialization_method == "binary_guided"
        population = binary_guided(environment)
    else
        println("Unknow initialization method ", environment.initialization_method)
    end
    return population
end