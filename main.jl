using Random
include("genetic_algorithm.jl")
#include("initialization.jl")
#include("individual.jl")
#include("selection.jl")
#include("mutation.jl")
#include("crossover.jl")

function main()
    println("Starting GA.")
    Random.seed!(0)

    A = rand(50,50) * ((10 - (0)) - 1)
    #println(A)

    ga = GeneticAlgorithm(1, "teste", 1000, 60, 30, 50, 50, 25, something, "binary", "random", "ranking", "fullyrandom", "binary_singlepoint", 0.1, "binary_mask", 0.8, 0.2, 0.5)
    #results = loop(ga)
    
    #println("Best fitness ", results[1].fitness)
    #println("Best solution ", results[1].chromosome)
end

main()