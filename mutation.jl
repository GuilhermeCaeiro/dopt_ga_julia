using StatsBase
using Distributions

function binary_singlepoint_mutation!(individual::Individual)
    ones = findall(x -> x == 1, individual.chromosome)
    zeros = findall(x -> x == 0, individual.chromosome)

    one_index = sample(ones, 1, replace=false)
    zero_index = sample(zeros, 1, replace=false)
    
    individual.chromosome[one_index[1]] = 0
    individual.chromosome[zero_index[1]] = 1
end

function binary_percentchange_mutation!(individual::Individual, total_ones::Int64, percent::Float64)
    ones = findall(x -> x == 1, individual.chromosome)
    zeros = findall(x -> x == 0, individual.chromosome)

    num_changes = ceil(total_ones * percent)
    ones_to_change = sample(ones, num_changes, replace=false)
    zeros_to_change = sample(zeros, num_changes, replace=false)

    individual.chromosome[ones_to_change] .= 0
    individual.chromosome[zeros_to_change] .= 1
end

function binary_variablepercentchange_mutation!(individual::Individual, total_ones::Int64, min_percent::Float64, max_percent::Float64)
    percent = rand(Uniform(min_percent, max_percent))
    binary_percentchange_mutation!(individual, total_ones, percent)
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