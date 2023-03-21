function LSFI(chromosome::Vector{Int64}, current_cost::Float64, Z::Matrix{Float64}, A::Matrix{Float64}, s::Int64)

    if isnothing(Z)
        try
            Z = inv(environment.A' * spdiagm(individual.chromosome) * environment.A)
        catch
            throw("Error: Z is nothing and cannot be calculated")
        end
    end

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