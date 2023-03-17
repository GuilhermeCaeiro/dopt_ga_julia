using Random
using Plots

include("utils.jl")
include("environment.jl")
include("individual.jl")
include("initialization.jl")
include("selection.jl")
include("mutation.jl")
include("crossover.jl")
# include("adaptation.jl")
include("genetic_algorithm.jl")

global current_generation = 0

function main()
    println("Starting GA.")
    Random.seed!(0)

    A, R, m, n, s = read_instance(300, 1)


    start_time = get_time_in_ms()

    println("Julia nthreads: ", Threads.nthreads())

    environment = Environment(
        1, # seed
        "teste", # instance
        10_000, # max_generations
        60, # max_time
        100, # population_size
        n, # n
        m, # m
        s, # s
        A, # A
        "binary", # encoding 
        "binary_random", # initialization_method
        [], # initialization_params
        "fullyrandom", # selecion_method
        "fullyrandom", # parent_selection_method
        [], # selection_params
        "binary_singlepoint", # mutation_method
        0.1, # mutation_probability
        [], # mutation_params
        "binary_singlepoint", # crossover_method
        0.9, # crossover_probability
        [], # crossover_params
        0.3, # elite_size
        0.5, # offspring_size
        "none", # adaptation_method | accepts "none" and "reset"
        [], # adaptation_params
        500, # generations_until_adaptation
        true, # perform path-relinking-like crossover
    )

    ga = GeneticAlgorithm(environment)
    # println(ga.environment.crossover_method)
    # println(ga.environment.mutation_method)
    # println(ga.environment.parent_selection_method)
    results, solutions = loop(ga)

    println("Num. generations: ", environment.max_generations, " Total Time: ", get_time_in_ms() - start_time)
    
    println("Fitness: ", results[1].fitness, " Objective Function: ", results[1].objective_function, " Penalty: ", results[1].penalty)
    println("Best chromosome, with ", sum(results[1].chromosome), " ones :", results[1].chromosome)
    println("Original fitness calculation: ", calculate_fitness(results[1].chromosome, environment.A, environment.s))

    savefig(plot(1:size(solutions, 1), solutions), "plot.png") 
    savefig(plot(1:length(ga.population_fitness_avg), ga.population_fitness_avg), "avg.png") 
    savefig(plot(1:length(ga.population_fitness_std), ga.population_fitness_std), "std.png") 
    savefig(plot(1:length(ga.population_fitness_avg_last_n_generations), ga.population_fitness_avg_last_n_generations), "avg_last_n.png") 
    savefig(plot(1:length(ga.population_fitness_std_last_n_generations), ga.population_fitness_std_last_n_generations), "avg_std_last_n.png")
    #println("Best fitness ", results[1].fitness)
    #println("Best solution ", results[1].chromosome)

    println("current_generation", current_generation)
end

main()