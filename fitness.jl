function calculate_objective_function(new_one_position::Int64, new_zero_position::Int64, current_cost::Float64, Z::Matrix{Float64}, A::Matrix{Float64}, det_tol::Float64=1.0e-10, newcost_tol::Float64=1.0e-6)
      
    w = A[new_one_position, :]
    r = A[new_zero_position, :]
    new_cost = -Inf
    det_sign = -1
    new_Z = nothing
    
    # calculating new cost
    wZ = w' * Z
    wZw = dot(wZ, w)
    rZr = dot(r, Z, r)
    wZr = dot(wZ, r)

    det_value = (1 + wZw) * (1 - rZr) + wZr^2
    
    if det_value > det_tol
        _new_cost = current_cost + log(det_value)
        if exp(_new_cost) > newcost_tol
            new_cost = _new_cost
            det_sign = 1
        end
    end

    # updating Z
    Y = Z - (wZ' * wZ) / (1 + dot(wZ, w))
    rY = r' * Y;
    new_Z = Y + (rY' * rY) / (1 - dot(rY, r))
    
    return new_cost, new_Z, det_sign
end

function ldet(A)
    (value,sign_det) = logabsdet(A);
    if (sign_det > 0) && (exp(value) > 1.0e-6)
        return value
    else
        return -Inf
    end
end

function calculate_fitness(chromosome, A, s)
    #println("CONVENTIONAL METHOD CALLED!")
    execution_statistics["conventional_of_calc_calls"][execution_statistics["current_generation"]] += 1

    num_ones = sum(chromosome)
    objective_function = ldet(A'*spdiagm(vec(chromosome))*A);
    penalty = - 100 * abs(num_ones - s)
    fitness = objective_function + penalty
    return fitness, objective_function, penalty
end

function calculate_fitness(chromosome::Vector{Int64}, A::Matrix{Float64}, s::Int64, parent::Individual, new_ones::Vector{Int64}, new_zeros::Vector{Int64})
    execution_statistics["efficient_of_calc_calls"][execution_statistics["current_generation"]] += 1
    
    Z = deepcopy(parent.Z_matrix)
    objective_function = deepcopy(parent.objective_function)
    det_sign = nothing
    
    actual_chromosome = deepcopy(parent.chromosome)
    extra_chromosome = nothing
    extra_fitness = -Inf
    extra_objfunc = deepcopy(parent.objective_function)
    extra_penalty = 0
    extra_Z = nothing

    pairs = Queue{Int64}()
    num_pairs = length(new_ones)
    count_iterations = 0

    for i in 1:num_pairs
        enqueue!(pairs, i)
    end

    while length(pairs) > 0
        count_iterations += 1
        i = dequeue!(pairs)
        actual_chromosome[new_ones[i]] = 1
        actual_chromosome[new_zeros[i]] = 0
        _objective_function, _Z, det_sign = calculate_objective_function(new_ones[i], new_zeros[i], objective_function, Z, A)

        if det_sign < 0
            #println("Sign < 0")
            if length(pairs) == 0
                #println("Final pair. The final OF ended up -inf.")
                objective_function = -Inf
                Z = nothing
            else
                if count_iterations <= num_pairs
                    #println("Pair ", new_ones[i], " ", new_zeros[i], " sent to de end of of the queue.")
                    enqueue!(pairs, i)
                    actual_chromosome[new_ones[i]] = 0
                    actual_chromosome[new_zeros[i]] = 1
                else
                    #println("Infinte objective function value. Falling back to the normal calculation method.")
                    fitness, objective_function, penalty = calculate_fitness(chromosome, A, s)
                    Z = nothing
                    break
                end
            end
        else
            objective_function = _objective_function
            Z = _Z

            if objective_function > extra_objfunc
                if parent.environment.perform_prlike_crossover
                    extra_Z = deepcopy(_Z)
                    extra_chromosome = deepcopy(actual_chromosome)
                    extra_objfunc = deepcopy(_objective_function)
                end
            end
        end
    end

    penalty = - 100 * abs(sum(chromosome) - s)
    fitness = objective_function + penalty

    if extra_objfunc > parent.objective_function
        extra_penalty = - 100 * abs(sum(extra_chromosome) - s)
        extra_fitness = extra_objfunc + extra_penalty
    else
        extra_fitness = -Inf
    end 

    return fitness, objective_function, penalty, Z, extra_chromosome, extra_fitness, extra_objfunc, extra_penalty, extra_Z
end
