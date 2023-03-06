using LinearAlgebra
using SparseArrays
using DataStructures

#include("genetic_algorithm.jl")

mutable struct Individual
    environment::Any
    chromosome::Vector{Int64}
    fitness::Float64
    objective_function::Float64
    penalty::Float64
    #parent::Individual
    #ones_changed::Vector{Int64}
    #zeros_changed::Vector{Int64}
    Z_matrix::Any
    
    function Individual(environment::Environment, chromosome::Vector{Int64})
        fitness, objective_function, penalty = calculate_fitness(chromosome, environment.A, environment.s)
        new(environment, chromosome, fitness, objective_function, penalty)
    end

    function Individual(environment::Environment, chromosome::Vector{Int64}, parent::Individual, new_ones::Vector{Int64}, new_zeros::Vector{Int64})
        if length(new_ones) != length(new_zeros)
            throw(error("Attempt to create Individual with length(ones_changed) != length(zeros_changed). This might indicate that the current solution (chromosome) is unfeasible."))
        end

        if !isdefined(parent, :Z_matrix) || isnothing(parent.Z_matrix)
            #println("UNDEFINED!!!!!!!!!!!!!!!!!!")
            try
                parent.Z_matrix = inv(environment.A' * diagm(parent.chromosome) * environment.A)

            catch e
                if isa(e, SingularException)
                    println("Matrix (environment.A' * diagm(parent.chromosome) * environment.A) is singular. Falling back to the normal calculation method.")
                    fitness, objective_function, penalty = calculate_fitness(chromosome, environment.A, environment.s)
                    return new(environment, chromosome, fitness, objective_function, penalty)
                else
                    throw(e)
                end
            end
        end

        if isinf(parent.objective_function)
            #println("Parent has infinete objective function value. Falling back to the normal calculation method.")
            fitness, objective_function, penalty = calculate_fitness(chromosome, environment.A, environment.s)
            new(environment, chromosome, fitness, objective_function, penalty)
        else
            fitness, objective_function, penalty, Z, extra_chromosome, extra_fitness, extra_objfunc, extra_penalty, extra_Z = calculate_fitness(
                chromosome, 
                environment.A, 
                environment.s,
                parent,
                new_ones,
                new_zeros
            )
            if isfinite(extra_fitness)
                return [new(environment, chromosome, fitness, objective_function, penalty, Z), new(environment, extra_chromosome, extra_fitness, extra_objfunc, extra_penalty, extra_Z)]
            else
                return new(environment, chromosome, fitness, objective_function, penalty, Z)
            end
        end

        
    end
end

function calculate_objective_function(new_one_position::Int64, new_zero_position::Int64, current_cost::Float64, Z, A)
    w = A[new_one_position, :]
    r = A[new_zero_position, :]
    new_cost = -Inf
    det_sign = -1

    # calculating new cost
    wZ = w' * Z
    wZw = dot(wZ, w)
    rZr = dot(r, Z, r)
    wZr = dot(wZ, r)

    det_value = (1 + wZw) * (1 - rZr) + wZr^2
    #det_value = (1 + wZ * w) * (1 - r' * Z * r) + (wZ * r)^2

    #println("det_value >>>>>>>>>>>>>>>>>>>>> ", det_value)
    
    #if det_value > 1.0e-6
    #    #println("new method's det > 0, ", det_value, ", log of that: ", log(det_value), " previous current cost", current_cost)
    #    new_cost = current_cost + log(det_value)
    #    det_sign = 1
    #end
    
    if det_value > 1.0e-10
        _new_cost = current_cost + log(det_value)
        if exp(_new_cost) > 1.0e-6
            new_cost = _new_cost
            det_sign = 1
        end
    end
    #    return current_cost, Z, -1
    #end

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
    num_ones = sum(chromosome)
    objective_function = ldet(A'*spdiagm(vec(chromosome))*A);
    #println("det normal fitness ", det(A'*spdiagm(vec(chromosome))*A))
    #println("log of det, normal fitness ", log(det(A'*spdiagm(vec(chromosome))*A)))
    #println(objective_function)
    penalty = - 100 * abs(num_ones - s)
    fitness = objective_function + penalty
    #println(chromosome)
    #println(num_ones, objective_function, penalty, fitness)

    """seeked_chromosome = [
        [1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1],
        [1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1],
        [1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1]
    ]

    if chromosome in seeked_chromosome
        
        println("Objective function calculated through the conventional fitness calculation. Value: ", objective_function, " Penalty: ", penalty)
        println(chromosome)
    end"""

    return fitness, objective_function, penalty
end

function calculate_fitness(chromosome::Vector{Int64}, A::Matrix{Float64}, s::Int64, parent::Individual, new_ones::Vector{Int64}, new_zeros::Vector{Int64})
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

    part_new_chromosome = deepcopy(parent.chromosome)
    """seeked_chromosome = [
        [1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1],
        [1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1],
        [1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1]
    ]"""

    #println((new_ones, chromosome[new_ones[1]], parent.chromosome[new_ones[1]], new_zeros, chromosome[new_zeros[1]], parent.chromosome[new_zeros[1]], num_pairs))

    for i in 1:num_pairs
        enqueue!(pairs, i)
    end

    while length(pairs) > 0
        count_iterations += 1
        i = dequeue!(pairs)
        actual_chromosome[new_ones[i]] = 1
        actual_chromosome[new_zeros[i]] = 0
        _objective_function, _Z, det_sign = calculate_objective_function(new_ones[i], new_zeros[i], objective_function, Z, A)
        #println("OF ", _objective_function, " sign ", det_sign)

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

        #if chromosome in seeked_chromosome
        #part_new_chromosome[new_ones[i]] = 1
        #part_new_chromosome[new_zeros[i]] = 0
        #println("pair (new one/new zero) ", new_ones[i], " ", new_zeros[i], " partial OF: ", _objective_function, " expected OF ", logabsdet(A'*spdiagm(vec(part_new_chromosome))*A)[1])

        #if objective_function != _objective_function
        #    part_new_chromosome[new_ones[i]] = 0
        #    part_new_chromosome[new_zeros[i]] = 1
        #end
        #end"""
    end

    reference_fitness, reference_of, reference_penalty = calculate_fitness(chromosome, A, s)
    
    if (!isapprox(reference_of, objective_function, atol = 0.00001)) #&& chromosome in seeked_chromosome
        println("Something went wrong! Result doesn't meet the reference.")
        println(" objective_function ", objective_function, " determinant of OF ", exp(objective_function), " reference_OF ", reference_of)
        println("chromosome ", chromosome)
        println("parent chromosome ", parent.chromosome)
        println("parent OF ", parent.objective_function)
        println("0-1 pairs: ", length(new_ones))
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
