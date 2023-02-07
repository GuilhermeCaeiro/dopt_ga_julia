using StatsBase

function binary_singlepoint_mutation!(individual::Individual)
    ones = findall(x -> x == 1, individual.chromosome)
    zeros = findall(x -> x == 0, individual.chromosome)

    one_index = sample(ones, 1, replace=false)
    zero_index = sample(zeros, 1, replace=false)
    
    individual.chromosome[one_index[1]] = 0
    individual.chromosome[zero_index[1]] = 1
end

function mutate(environment::Environment, individual::Individual)
    new_individual = deepcopy(individual)
    
    if environment.mutation_method == "binary_singlepoint"
        binary_singlepoint_mutation!(new_individual)
    elseif environment.mutation_method == "something"
        #
    else
        println("Unknow mutation_method method ", environment.mutation_method)
    end
    return new_individual
end