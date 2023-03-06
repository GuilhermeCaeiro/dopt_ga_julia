using StatsBase

#include("genetic_algorithm.jl")

function binary_singlepoint(chromosome1::Vector{Int64}, chromosome2::Vector{Int64}, num_ones::Integer)
    cut_position = sample(1:size(chromosome1, 1), 1, replace=false)[1]
    #println((cut_position, typeof(cut_position)))
    #println((chromosome1, chromosome2))
    #cut_position = 1 #size(chromosome1, 1)


    if cut_position == 1 || cut_position == size(chromosome1, 1)
        #println("cut_position == 0 || cut_position == size(chromosome1, 1)")
        return chromosome1, chromosome2
    else
        offspring1 = [chromosome1[1:cut_position]; chromosome2[(cut_position + 1):end]]
        offspring2 = [chromosome2[1:cut_position]; chromosome1[(cut_position + 1):end]]
        
        offspring1_num_ones = sum(offspring1)
        offspring2_num_ones = sum(offspring2)

        #println((num_ones, size(offspring1, 1), offspring1_num_ones, size(offspring2, 1), offspring2_num_ones, offspring1, offspring2))

        if offspring1_num_ones < num_ones
            #println("offspring1_num_ones < num_ones")
            zeros = findall(x -> x == 0, offspring1)
            selected_positions = sample(zeros, num_ones - offspring1_num_ones, replace=false)
            #println(selected_positions)
            offspring1[selected_positions] .= 1
        elseif offspring1_num_ones > num_ones
            #println("offspring1_num_ones > num_ones")
            ones = findall(x -> x == 1, offspring1)
            selected_positions = sample(ones, offspring1_num_ones - num_ones, replace=false)
            #println(selected_positions)
            offspring1[selected_positions] .= 0
        end

        if offspring2_num_ones < num_ones
            #println("offspring2_num_ones < num_ones")
            zeros = findall(x -> x == 0, offspring2)
            #println("zeros", zeros)
            selected_positions = sample(zeros, num_ones - offspring2_num_ones, replace=false)
            #println(selected_positions)
            offspring2[selected_positions] .= 1
        elseif offspring2_num_ones > num_ones
            #println("offspring2_num_ones > num_ones")
            ones = findall(x -> x == 1, offspring2)
            selected_positions = sample(ones, offspring2_num_ones - num_ones, replace=false)
            #println(selected_positions)
            offspring2[selected_positions] .= 0
        end

        #println((size(offspring1, 1), sum(offspring1), size(offspring2, 1), sum(offspring2), offspring1, offspring2))

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
    #offspring = []
    p1_chromosome = deepcopy(parent1.chromosome) 
    p2_chromosome = deepcopy(parent2.chromosome)
    chromosome1 = []
    chromosome2 = []

    if environment.crossover_method == "binary_singlepoint"
        chromosome1, chromosome2 = binary_singlepoint(p1_chromosome, p2_chromosome, environment.s)
        #offspring = [Individual(environment, chromosome1), Individual(environment, chromosome2)]
    elseif environment.crossover_method == "binary_mask"
        #println("type of chromosome: ", typeof(parent1.chromosome))
        chromosome1, chromosome2 = binary_mask(p1_chromosome, p2_chromosome, environment.s)
        #offspring = [Individual(environment, chromosome1), Individual(environment, chromosome2)]
    else
        throw(error("Uknown crossover method ", environment.crossover_method))
    end

    differ_par1_offs1 = p1_chromosome - chromosome1
    differ_par1_offs2 = p1_chromosome - chromosome2
    differ_par2_offs1 = p2_chromosome - chromosome1
    differ_par2_offs2 = p2_chromosome - chromosome2

    val_differ_par1_offs1 = sum(abs.(differ_par1_offs1))
    val_differ_par1_offs2 = sum(abs.(differ_par1_offs2))
    val_differ_par2_offs1 = sum(abs.(differ_par2_offs1))
    val_differ_par2_offs2 = sum(abs.(differ_par2_offs2))

    parent_offs1 = parent1
    parent_offs2 = parent2
    
    new_ones_offs1 = []
    new_ones_offs2 = []
    new_zeros_offs1 = []
    new_zeros_offs2 = []

    if (val_differ_par1_offs1 >= val_differ_par2_offs1)
        parent_offs1 = parent2
        new_ones_offs1 = findall(x -> x == -1, differ_par2_offs1)
        new_zeros_offs1 = findall(x -> x == 1, differ_par2_offs1)
    else
        parent_offs1 = parent1
        new_ones_offs1 = findall(x -> x == -1, differ_par1_offs1)
        new_zeros_offs1 = findall(x -> x == 1, differ_par1_offs1)
    end

    if (val_differ_par1_offs2 >= val_differ_par2_offs2)
        parent_offs2 = parent2
        new_ones_offs2 = findall(x -> x == -1, differ_par2_offs2)
        new_zeros_offs2 = findall(x -> x == 1, differ_par2_offs2)
    else
        parent_offs2 = parent1
        new_ones_offs2 = findall(x -> x == -1, differ_par1_offs2)
        new_zeros_offs2 = findall(x -> x == 1, differ_par1_offs2)
    end

    

    offspring = Vector{Individual}()
    ind1 = Individual(environment, chromosome1, parent_offs1, new_ones_offs1, new_zeros_offs1)
    ind2 = Individual(environment, chromosome2, parent_offs2, new_ones_offs2, new_zeros_offs2)
    offspring = [offspring; ind1; ind2]

    return offspring
end

#println("Hello")