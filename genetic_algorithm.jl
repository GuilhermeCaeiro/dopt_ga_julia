using Random

#include("individual.jl")

mutable struct GeneticAlgorithm
    environment::Environment
    elite_size::Int64
    offspring_size::Int64
    best_solution::Float64
    best_solution_tracking::Vector{Float64}
    elite::Vector{Individual}

    function GeneticAlgorithm(environment::Environment)
        #set seed
        Random.seed!(environment.seed)

        elitesize = ceil(Int64, environment.population_size * environment.elite_size)
        offspringsize = floor(Int64, environment.population_size * environment.offspring_size)
        bestsolution = -Inf
        elite = Vector{Individual}()
        best_solution_tracking = Vector{Float64}()

        # populate struct variables
        new(
            environment,
            elitesize,
            offspringsize,
            bestsolution,
            best_solution_tracking,
            elite
        )

    end

end

function loop(ga::GeneticAlgorithm)
    # Initialize population
    #population = Vector{Individual}()
    population = initialize_population(ga.environment)

    # runs generations
    for generation in 1:ga.environment.max_generations
        if generation % 1000 == 0
            println("Generation ", generation)
        end
        sort!(population, by = v -> v.fitness, rev = true)
        
        children = Vector{Individual}()
        
        # crossover and mutation
        while size(children, 1) < ga.offspring_size
            if rand(Float64, 1)[1] < ga.environment.crossover_probability
                parents = select(ga.environment, population, 2, ga.environment.parent_selection_method)
                #print(methods(breed))
                offspring = breed(ga.environment, parents[1], parents[2])
                children = [children; offspring]
            end

            if rand(Float64, 1)[1] < ga.environment.mutation_probability
                somebody = select(ga.environment, population, 1, ga.environment.parent_selection_method)
                mutant = mutate(ga.environment, somebody[1])
                children = [children; mutant]
            end
        end

        population = [population; children]
        sort!(population, by = v -> v.fitness, rev = true)

        if population[1].fitness > ga.best_solution
            ga.best_solution = population[1].fitness
        end

        push!(ga.best_solution_tracking, ga.best_solution)

        #println("ga.elite_size ", ga.elite_size)
        ga.elite = population[1:ga.elite_size]
        commoners = population[(ga.elite_size + 1):end]
        commoners = select(ga.environment, commoners, ga.environment.population_size - size(ga.elite)[1], ga.environment.selecion_method)
        population = [ga.elite; commoners]
    end

    return ga.elite, ga.best_solution_tracking
end