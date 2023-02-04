using StatsBase

function binary_singlepoint_mutation!(individual::Individual)
    ones = findall(x -> x == 1, v)
    zeros = findall(x -> x == 0, v)

    one_index = sample(ones, 1, replace=false)
    zero_index = sample(zeros, 1, replace=false)
    
    individual[one_index] = 0
    individual[zero_index] = 1
end

function mutate(ga::GeneticAlgorithm, individual::Individual)
    new_individual = deepcopy(individual)
    
    if ga.mutation_method == "binary_singlepoint"
        binary_singlepoint_mutation!(new_individual)
    elseif ga.mutation_method == "something"
        #
    else
        println("Unknow mutation_method method ", ga.mutation_method)
    end
    return new_individual
end