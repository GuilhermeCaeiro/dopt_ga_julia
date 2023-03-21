using Random
using Dates
using Statistics

#include("individual.jl")

mutable struct GeneticAlgorithm
    environment::Environment
    elite_size::Int64
    offspring_size::Int64
    best_solution::Float64
    best_solution_tracking::Vector{Float64}
    best_solution_change_or_adaptation_performed::Int64
    elite::Vector{Individual}
    population_fitness_avg::Vector{Float64}
    population_fitness_std::Vector{Float64}
    population_fitness_avg_last_n_generations::Vector{Float64}
    population_fitness_std_last_n_generations::Vector{Float64}

    function GeneticAlgorithm(environment::Environment)
        #set seed
        Random.seed!(environment.seed)

        elitesize = ceil(Int64, environment.population_size * environment.elite_size)
        offspringsize = floor(Int64, environment.population_size * environment.offspring_size)
        bestsolution = -Inf
        elite = Vector{Individual}()
        best_solution_tracking = Vector{Float64}()
        population_fitness_avg = Vector{Float64}()
        population_fitness_std = Vector{Float64}()
        population_fitness_avg_last_n_generations = Vector{Float64}()
        population_fitness_std_last_n_generations = Vector{Float64}()

        # populate struct variables
        new(
            environment,
            elitesize,
            offspringsize,
            bestsolution,
            best_solution_tracking,
            1,
            elite,
            population_fitness_avg,
            population_fitness_std,
            population_fitness_avg_last_n_generations,
            population_fitness_std_last_n_generations
        )

    end

end

function loop(ga::GeneticAlgorithm)
    # Initialize population
    population = initialize_population(ga.environment)
    #iter_times = Vector{Int64}()

    println("Post initialization stats: Conv. ->", execution_statistics["conventional_of_calc_calls"], " Effic. ->", execution_statistics["efficient_of_calc_calls"])

    loop_start_time = get_time_in_ms()
    # runs generations
    for generation in 1:ga.environment.max_generations
        execution_statistics["current_generation"] = generation
        execution_statistics["conventional_of_calc_calls"][generation] = 0
        execution_statistics["efficient_of_calc_calls"][generation] = 0
        execution_statistics["parent_zcalc_from_children"][generation] = 0
        iter_start_time = get_time_in_ms()

        sort!(population, by = v -> v.fitness, rev = true)

        children = Vector{Individual}()

        # crossover and mutation
        while size(children, 1) < ga.offspring_size
            if rand(Float64, 1)[1] < ga.environment.crossover_probability
                parents = select(ga.environment, population, 2, ga.environment.parent_selection_method)
                offspring = breed(ga.environment, parents[1], parents[2])
                children = [children; offspring]
            end

            if rand(Float64, 1)[1] < ga.environment.mutation_probability
                somebody = select(ga.environment, population, 1, ga.environment.parent_selection_method)
                mutant = mutate(ga.environment, somebody[1])
                children = [children; mutant]
            end
        end

        # ======================================================
        # SUGESTTION FOR ADAPTATION
        # old_mutation_method = ga.environment.mutation_method
        # ga.environment.mutation_method = "binary_search"
        # mutant = mutate(ga.environment, population[1])
        # children = [children; mutant]
        # ga.environment.mutation_method = old_mutation_method
        # ======================================================

        population = [population; children]
        sort!(population, by = v -> v.fitness, rev = true)

        if population[1].fitness > ga.best_solution
            ga.best_solution = population[1].fitness
            ga.best_solution_change_or_adaptation_performed = generation
        end

        push!(ga.best_solution_tracking, ga.best_solution)

        ga.elite = population[1:ga.elite_size]
        commoners = population[(ga.elite_size + 1):end]
        commoners = select(ga.environment, commoners, ga.environment.population_size - ga.elite_size, ga.environment.selection_method)
        population = [ga.elite; commoners]

        println(
            "Iteration: ", generation, 
            " Pop. size: ", length(population), 
            " Best Sol.: ", population[1].objective_function, 
            " Avg Sol.: ", mean(x-> x.fitness, population), 
            " Std Sol.: ", std([x.fitness for x in population]), 
            " Conventional OF: ", execution_statistics["conventional_of_calc_calls"][generation], 
            " Efficient OF: ", execution_statistics["efficient_of_calc_calls"][generation],
            " Par. Zcalc f/ children: ", execution_statistics["parent_zcalc_from_children"][generation]
        )

        iter_time = get_time_in_ms() - iter_start_time

        push!(execution_statistics["iter_times"], iter_time.value)

        population_fitness = [individual.fitness for individual in population]
        push!(ga.population_fitness_avg, mean(population_fitness))
        push!(ga.population_fitness_std, std(population_fitness))


        start_position = 1
        if generation > ga.environment.generations_until_adaptation
            start_position = generation - ga.environment.generations_until_adaptation
        end
        
        push!(ga.population_fitness_avg_last_n_generations, mean(ga.population_fitness_avg[start_position:end]))
        push!(ga.population_fitness_std_last_n_generations, std(ga.population_fitness_std[start_position:end]))
        
        if (((get_time_in_ms() - loop_start_time).value / 1000.0) >= ga.environment.max_time) # "/ 1000" because the time provided in the setup file should be in seconds.
            println("Maximum time reached. Stopping the GA.")
            break
        end

        if ga.environment.adaptation_method != "none"
            # checks if it is time to adapt
            if (generation - ga.best_solution_change_or_adaptation_performed) > ga.environment.generations_until_adaptation
                println("Performing adaptation in generation ", generation) 
                population = adapt(ga.environment, ga.elite)
                ga.best_solution_change_or_adaptation_performed = generation
            end
        end
    end

    println(
        "Avg. Time for ", ga.environment.max_generations, 
        " iterations: ", mean(execution_statistics["iter_times"]), 
        " Min|Max: ", minimum(execution_statistics["iter_times"]), " | ", maximum(execution_statistics["iter_times"])
    )

    return ga.elite, ga.best_solution_tracking
end