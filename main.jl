using Random
using Dates
using MAT
using DelimitedFiles
using Plots

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

function read_instance(m, i)
    file = matopen("instances/instance_$(m)_$(i).mat")
    A = read(file, "A")
    R = read(file, "R")
    close(file)
    n, m = size(A, 1), size(A, 2) # Gabriel uses n as the number of lines (experiments).
    s = Int(n/2)
    return A, R, m, n, s
end

function main()
    println("Starting GA.")
    Random.seed!(0)
    
    # loading instance
    #A = rand(50,50) * ((10 - (0)) - 1)
    #println(A)

    A, R, m, n, s = read_instance(40, 1)


    start_time = get_time_in_ms()

    environment = Environment(
        1, 
        "teste", 
        20000, 
        60, 
        30, 
        n, 
        m, 
        s, 
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
    results, solutions = loop(ga)

    println("Total Time: ", get_time_in_ms() - start_time)
    
    print(results[1].fitness)

    p = plot(1:size(solutions, 1), solutions)
    savefig(p, "plot.png") 
    #println("Best fitness ", results[1].fitness)
    #println("Best solution ", results[1].chromosome)
end

main()