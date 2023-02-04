using StatsBase

function binary_mask(chromosome1::Vector{Integer}, chromosome2::Vector{Integer}, num_ones::Integer)
    ones = unique([findall(x -> x == 1, chromosome1); findall(x -> x == 1, chromosome2)])
    chromosome_size = size(chromosome1, 1)

    offspring1 = zeros(chromosome_size)
    offspring2 = zeros(chromosome_size)

    offspring1[sample(ones, num_ones, replace=false)] .= 1
    offspring2[sample(ones, num_ones, replace=false)] .= 1

    return offspring1, offspring2

end

function breed(ga::GeneticAlgorith, parent1::Individual, parent2::Individual)
    offspring = []

    if ga.crossover_method == "binary_mask"
        chromosome1, chromosome2 = binary_mask(parent1.chromosome, parent2.chromosome, ga.s)
        offspring = [Individual(ga, chromosome1), Individual(ga, chromosome2)]
    else
        println("Uknown crossover method ", ga.crossover_method)
    end
    
    return offspring
end