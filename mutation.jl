using StatsBase
using Distributions

function binary_search_mutation(chromosome::Vector{Int64}, current_cost::Float64, Z::Matrix{Float64}, A::Matrix{Float64}, s::Int64)
    n = length(chromosome)
    ones = findall(x -> x == 1, chromosome)
    zeros = findall(x -> x == 0, chromosome)
    n2 = length(ones)
    new_one = -1
    new_zero = -1
    actual_cost = current_cost
    stop = false
    for i in 1:n2
        for j in i:n2
            new_cost, new_Z, det_sign = calculate_objective_function(zeros[i], ones[j], current_cost, Z, A)
            if det_sign > -1
                if new_cost > actual_cost
                    actual_cost = new_cost
                    new_one = zeros[i]
                    new_zero = ones[j]
                    Z = new_Z
                    stop = true
                    break
                end
            end
        end
        if stop
            break
        end
    end
    if new_one != -1
        chromosome[new_zero] = 0
        chromosome[new_one] = 1
    else
        one_index = sample(ones, 1, replace=false)
        zero_index = sample(zeros, 1, replace=false)
        chromosome[one_index[1]] = 0
        chromosome[zero_index[1]] = 1
        actual_cost, Z, det_sign = calculate_objective_function(zero_index[1], one_index[1], current_cost, Z, A)
    end
    # _fitness, _objective_function, _penalty = calculate_fitness(chromosome, A, s)
    # println("correct: ", _fitness, " calculated: ", actual_cost)
    return chromosome, actual_cost, Z
end

function binary_singlepoint_mutation(chromosome::Vector{Int64}, environment::Environment)
    ones = findall(x -> x == 1, chromosome)
    zeros = findall(x -> x == 0, chromosome)

    one_index = sample(ones, 1, replace=false)
    zero_index = sample(zeros, 1, replace=false)
    
    chromosome[one_index[1]] = 0
    chromosome[zero_index[1]] = 1

    # zero_index => position that became one
    # one_index => position that became zero
    return chromosome, zero_index, one_index 
end

function binary_percentchange_mutation(chromosome::Vector{Int64}, environment::Environment, percent)#total_ones::Int64, percent::Float64)
    ones = findall(x -> x == 1, chromosome)
    zeros = findall(x -> x == 0, chromosome)

    num_changes = ceil(Int64, environment.s * percent)
    ones_to_change = sample(ones, num_changes, replace=false)
    zeros_to_change = sample(zeros, num_changes, replace=false)

    chromosome[ones_to_change] .= 0
    chromosome[zeros_to_change] .= 1

    # zeros_to_change => positions that became one
    # ones_to_change => positions that became zero
    return chromosome, zeros_to_change, ones_to_change
end

function binary_variablepercentchange_mutation(chromosome::Vector{Int64}, environment::Environment, min_percent, max_percent) #total_ones::Int64, min_percent::Float64, max_percent::Float64)
    percent = rand(Uniform(min_percent, max_percent))
    return binary_percentchange_mutation(chromosome, environment, percent)
end

function mutate(environment::Environment, individual::Individual)
    new_chromosome = deepcopy(individual.chromosome)
    fitness = deepcopy(individual.fitness)

    if !isdefined(individual, :Z_matrix) || isnothing(individual.Z_matrix)
        Z = inv(environment.A' * spdiagm(individual.chromosome) * environment.A)
        individual.Z_matrix = deepcopy(Z)
    else
        Z = deepcopy(individual.Z_matrix)
    end
    new_ones = []
    new_zeros = []
    
    if environment.mutation_method == "binary_singlepoint"
        new_chromosome, new_ones, new_zeros = binary_singlepoint_mutation(new_chromosome, environment)
    elseif environment.mutation_method == "binary_percentchange"
        # mutation_params is assumed to be a float between 0 and 1.
        new_chromosome, new_ones, new_zeros = binary_percentchange_mutation!(new_chromosome, environment, environment.mutation_params[1])
    elseif environment.mutation_method == "binary_variablepercentchange"
        # mutation_params is assumed to be two floats between 0 and 1, with the fist smaller than the second.
        new_chromosome, new_ones, new_zeros = binary_variablepercentchange_mutation!(new_chromosome, environment, environment.mutation_params[1], environment.mutation_params[2])
    elseif environment.mutation_method == "binary_search"
        new_chromosome, new_fitness, new_Z = binary_search_mutation(new_chromosome, fitness, Z, environment.A, environment.s)
        return Individual(environment, new_chromosome, new_fitness, new_fitness, 0.0, new_Z)
    else
        throw(error("Unknow mutation_method method " + environment.mutation_method))
    end
    
    return Individual(environment, new_chromosome, individual, new_ones, new_zeros)
end