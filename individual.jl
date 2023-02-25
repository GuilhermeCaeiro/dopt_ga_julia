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
            parent.Z_matrix = inv(environment.A' * diagm(parent.chromosome) * environment.A)
        end

        if isinf(parent.objective_function)
            println("Parent has infinete objective function value. Falling back to the normal calculation method.")
            fitness, objective_function, penalty = calculate_fitness(chromosome, environment.A, environment.s)
            new(environment, chromosome, fitness, objective_function, penalty)
        else
            fitness, objective_function, penalty, Z = calculate_fitness(
                chromosome, 
                environment.A, 
                environment.s,
                parent,
                new_ones,
                new_zeros
            )

            new(environment, chromosome, fitness, objective_function, penalty, Z)
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

    println("det_value ", det_value)
    
    if det_value > 0
        println("new method's det > 0, ", det_value, ", log of that: ", log(det_value), " previous current cost", current_cost)
        new_cost = current_cost + log(det_value)
        det_sign = 1
    end
    #    return current_cost, Z, -1
    #end

    # updating Z
    Y = Z - (wZ' * wZ) / (1 + wZ * w)
    rY = r' * Y;
    new_Z = Y + (rY' * rY) / (1 - rY * r)
    
    return new_cost, new_Z, det_sign
end

function ldet(A)
    (value,sign_det) = logabsdet(A);
    if (sign_det > 0)
        return value
    else
        return -Inf
    end
end

function calculate_fitness(chromosome, A, s)
    num_ones = sum(chromosome)
    objective_function = ldet(A'*spdiagm(vec(chromosome))*A);
    println("det normal fitness ", det(A'*spdiagm(vec(chromosome))*A))
    #println("log of det, normal fitness ", log(det(A'*spdiagm(vec(chromosome))*A)))
    println(objective_function)
    penalty = - 100 * abs(num_ones - s)
    fitness = objective_function + penalty
    #println(chromosome)
    #println(num_ones, objective_function, penalty, fitness)
    return fitness, objective_function, penalty
end

function calculate_fitness(chromosome::Vector{Int64}, A::Matrix{Float64}, s::Int64, parent::Individual, new_ones::Vector{Int64}, new_zeros::Vector{Int64})
    Z = deepcopy(parent.Z_matrix)
    objective_function = deepcopy(parent.objective_function)

    pairs = Queue{Int64}()
    num_pairs = length(new_ones)
    count_iterations = 0

    println((new_ones, chromosome[new_ones[1]], parent.chromosome[new_ones[1]], new_zeros, chromosome[new_zeros[1]], parent.chromosome[new_zeros[1]], num_pairs))

    for i in 1:num_pairs
        enqueue!(pairs, i)
    end

    while length(pairs) > 0
        count_iterations += 1
        i = dequeue!(pairs)
        _objective_function, _Z, det_sign = calculate_objective_function(new_ones[i], new_zeros[i], objective_function, Z, A)
        println("OF ", _objective_function, " sign ", det_sign)

        if det_sign < 0
            println("Sign < 0")
            if length(pairs) == 0
                objective_function = -Inf
                Z = nothing
            else
                if count_iterations <= num_pairs
                    enqueue!(pairs, i)
                else
                    println("Infinte objective function value. Falling back to the normal calculation method.")
                    fitness, objective_function, penalty = calculate_fitness(chromosome, A, s)
                    return fitness, objective_function, penalty, nothing
                end
            end
        else
            objective_function = _objective_function
            Z = _Z
        end
    end

    num_ones = sum(chromosome)
    penalty = - 100 * abs(num_ones - s)
    fitness = objective_function + penalty

    return fitness, objective_function, penalty, Z
end
