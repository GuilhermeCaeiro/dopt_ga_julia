using Random
include("individual.jl")

mutable struct GeneticAlgorithm
    seed::Int64
    instance::String
    max_generations::Int64
    max_time::Int64
    population_size::Int64
    n::Int64
    m::Int64
    s::Int64
    A::Matrix{Float64}
    encoding::String
    initialization_method::String
    selecion_method::String
    parent_selection_method::String
    mutation_method::String
    mutation_probability::Float64
    crossover_method::String
    crossover_probability::Float64
    best_solution::Float64
    elite_size::Int64
    offspring_size::Int64
    elite::Vector{Individual}

    function GeneticAlgorithm(
        seed::Int64,
        instance::String,
        max_generations::Int64,
        max_time::Int64,
        population_size::Int64,
        n::Int64,
        m::Int64,
        s::Int64,
        A::Matrix{Float64},
        encoding::String,
        initialization_method::String,
        selecion_method::String,
        parent_selection_method::String,
        mutation_method::String,
        mutation_probability::Float64,
        crossover_method::String,
        crossover_probability::Float64,
        elite_size::Int64,
        offspring_size::Float64
    )
        #set seed
        Random.seed!(seed)

        elitesize = ceil(population_size * elite_size)
        offspringsize = floor(population_size * offspring_size)
        bestsolution = -Inf
        elite = Vector{Individual}()

        # populate struct variables
        new(
            seed,
            instance,
            max_generations,
            max_time,
            population_size,
            n,
            m,
            s,
            A,
            encoding,
            initialization_method,
            selecion_method,
            parent_selection_method,
            mutation_method,
            mutation_probability,
            crossover_method,
            crossover_probability,
            bestsolution,
            elitesize,
            offspringsize,
            elite
        )

    end

end

function loop(ga::GeneticAlgorithm)
    # Initialize population
    population = Vector{Individual}()

    

    # runs generations
    for generation in 1:ga.max_generations
        sort!(population, by = v -> v.fitness, rev = true)
        
        children = Vector{Individual}()
        
        # crossover and mutation
        while size(children, 1) < ga.offspring_size
            if rand(Float64, 1) < ga.crossover_probability
                parents = select(ga, population, 2, ga.parent_selection_method)
                offspring = breed(parents[1], parents[2])
                children = [children; offspring]
            end

            if rand(Float64, 1) < ga.mutation_probability
                somebody = select(ga, population, 1, ga.parent_selection_method)
                mutant = mutate(ga, somebody[1])
                children = [children; mutant]
            end
        end

        population = [population; children]
        sort!(population, by = v -> v.fitness, rev = true)

        if population[1].fitness > ga.best_solution
            ga.best_solution = population[1].fitness
        end

        ga.elite = population[1:ga.elite_size]
        commoners = population[(ga.elite_size + 1):end]
        commoners = select(ga, commoners, ga.population_size - size(ga.elite)[1], ga.selecion_method)
        population = [ga.elite; commoners]
    end

    return ga.elite
end