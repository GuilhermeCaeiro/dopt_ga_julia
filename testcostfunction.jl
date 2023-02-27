using SparseArrays
using LinearAlgebra
using Random
using StatsBase

include("utils.jl")
include("environment.jl")
include("individual.jl")

A, R, m, n, s = read_instance(300, 1)

environment = Environment(
    1, 
    "teste", 
    6, 
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

for j in 12:100
    Random.seed!(j)
    chromosome = zeros(Int64, environment.n)
    indices = sample(1:environment.n, environment.s, replace=false)
    chromosome[indices] .= 1

    individual = Individual(environment, chromosome)
    println("Individual OF: ", individual.objective_function)

    individual.Z_matrix = inv(environment.A' * spdiagm(chromosome) * environment.A)

    _ones = findall(x -> x == 1, chromosome)
    _zeros = findall(x -> x == 0, chromosome)

    

    for i in 1:length(_ones)
        new_chromosome = deepcopy(chromosome)
        new_chromosome[_ones[i]] = 0
        new_chromosome[_zeros[i]] = 1

        test_individual = Individual(environment, new_chromosome, individual, [_zeros[i]], [_ones[i]])
        control_individual = Individual(environment, new_chromosome)
        println("Parent OF: ", individual.objective_function, ", test child OF - parent's: ", test_individual.objective_function - individual.objective_function)
        println("Parent OF: ", individual.objective_function, ", control child OF - parent's: ", control_individual.objective_function - individual.objective_function)

        println((
            j,
            i,
            "Test", test_individual.objective_function, 
            "Control", control_individual.objective_function, 
            "Same", isapprox(test_individual.objective_function, control_individual.objective_function, atol=0.000001)
        ))

        if !isapprox(test_individual.objective_function, control_individual.objective_function, atol=0.000001)
            println((
            j,
            i,
            "Test", test_individual.objective_function, 
            "Control", control_individual.objective_function, 
            "Same", isapprox(test_individual.objective_function, control_individual.objective_function, atol=0.000001)
        ))
            throw(error("Different OF values found!"))
        end
    end
end