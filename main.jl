using Random
using Dates

include("environment.jl")
include("individual.jl")
include("initialization.jl")
include("selection.jl")
include("mutation.jl")
include("crossover.jl")
include("genetic_algorithm.jl")


function get_time_in_ms()
    return convert(Dates.Millisecond, Dates.now())
end

function main()
    println("Starting GA.")
    Random.seed!(0)

    start_time = get_time_in_ms()

    A = rand(50,50) * ((10 - (0)) - 1)
    #println(A)

    environment = Environment(
        1, 
        "teste", 
        1000, 
        60, 
        30, 
        50, 
        50, 
        25, 
        A, 
        "binary", 
        "binary_random", 
        "ranking", 
        "fullyrandom", 
        "binary_singlepoint", 
        0.1, 
        "binary_mask", 
        0.8, 
        0.2, 
        0.5
    )

    ga = GeneticAlgorithm(environment)
    println(ga.environment.crossover_method)
    results = loop(ga)

    println("Total Time: ", get_time_in_ms() - start_time)
    
    print(results)
    #println("Best fitness ", results[1].fitness)
    #println("Best solution ", results[1].chromosome)
end

main()