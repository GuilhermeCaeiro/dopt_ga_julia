using StatsBase

#include("genetic_algorithm.jl")

function binary_mask(chromosome1::Vector{Int64}, chromosome2::Vector{Int64}, num_ones::Integer)
    ones = unique([findall(x -> x == 1, chromosome1); findall(x -> x == 1, chromosome2)])
    chromosome_size = size(chromosome1, 1)

    offspring1 = zeros(Int64, chromosome_size)
    offspring2 = zeros(Int64, chromosome_size)

    offspring1[sample(ones, num_ones, replace=false)] .= 1
    offspring2[sample(ones, num_ones, replace=false)] .= 1

    #println("offspring1 ", offspring1)

    return offspring1, offspring2

end

function breed(environment::Environment, parent1::Individual, parent2::Individual)
    offspring = []

    if environment.crossover_method == "binary_mask"
        #println("type of chromosome: ", typeof(parent1.chromosome))
        chromosome1, chromosome2 = binary_mask(parent1.chromosome, parent2.chromosome, environment.s)
        offspring = [Individual(environment, chromosome1), Individual(environment, chromosome2)]
    else
        println("Uknown crossover method ", environment.crossover_method)
    end
    
    return offspring
end

#println("Hello")