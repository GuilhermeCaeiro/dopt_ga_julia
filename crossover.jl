using StatsBase

#include("genetic_algorithm.jl")

function binary_singlepoint(chromosome1::Vector{Int64}, chromosome2::Vector{Int64}, num_ones::Integer)
    cut_position = sample(1:size(chromosome1, 1), 1, replace=false)

    if cut_position == 0
        return deepcopy(chromosome1), deepcopy(chromosome2)
    else
        offspring1 = [chromosome1[1:cut_position]; chromosome2[(cut_position + 1):end]]
        offspring1 = [chromosome2[1:cut_position]; chromosome1[(cut_position + 1):end]]
        
        offspring1_num_ones = sum(offspring1)
        offspring2_num_ones = sum(offspring2)

        if offspring1_num_ones < num_ones
            zeros = findall(x -> x == 0, chromosome1)
            selected_positions = sample(zeros, num_ones - offspring1_num_ones, replace=false)
            offspring1[selected_positions] .= 1
        elseif offspring1_num_ones > num_ones
            ones = findall(x -> x == 1, chromosome1)
            selected_positions = sample(ones, offspring1_num_ones - num_ones, replace=false)
            offspring1[selected_positions] .= 0
        end

        if offspring2_num_ones < num_ones
            zeros = findall(x -> x == 0, chromosome2)
            selected_positions = sample(zeros, num_ones - offspring2_num_ones, replace=false)
            offspring2[selected_positions] .= 1
        elseif offspring2_num_ones > num_ones
            ones = findall(x -> x == 1, chromosome2)
            selected_positions = sample(ones, offspring2_num_ones - num_ones, replace=false)
            offspring2[selected_positions] .= 0
        end

        return offspring1, offspring2
    end
end

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