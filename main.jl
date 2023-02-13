using Random
using Plots

include("utils.jl")
include("environment.jl")
include("individual.jl")
include("initialization.jl")
include("selection.jl")
include("mutation.jl")
include("crossover.jl")
include("genetic_algorithm.jl")

function main()
    println("Starting GA.")
    Random.seed!(0)

    A, R, m, n, s = read_instance(40, 1)


    start_time = get_time_in_ms()

    environment = Environment(
        1, 
        "teste", 
        10000, 
        60, 
        100, 
        n, 
        m, 
        s, 
        A, 
        "binary", 
        "binary_random", 
        [],
        "roulette", 
        "fullyrandom", 
        [],
        "binary_singlepoint", 
        0.1, 
        [],
        "binary_mask", 
        0.9, 
        [],
        0.3, 
        0.5
    )

    ga = GeneticAlgorithm(environment)
    println(ga.environment.crossover_method)
    results, solutions = loop(ga)

    println("Num. generations: ", environment.max_generations, " Total Time: ", get_time_in_ms() - start_time)
    
    print(results[1].fitness)

    p = plot(1:size(solutions, 1), solutions)
    savefig(p, "plot.png") 
    #println("Best fitness ", results[1].fitness)
    #println("Best solution ", results[1].chromosome)
end

main()