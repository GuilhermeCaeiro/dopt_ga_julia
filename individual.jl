using LinearAlgebra
using SparseArrays
using DataStructures

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
        #println("CONVENTIONAL METHOD CALLED!")
        execution_statistics["conventional_of_calc_calls"][execution_statistics["current_generation"]] += 1
        fitness, objective_function, penalty = calculate_fitness(chromosome, environment.A, environment.s)
        new(environment, chromosome, fitness, objective_function, penalty)
    end

    function Individual(environment::Environment, chromosome::Vector{Int64}, fitness::Float64, objective_function::Float64, penalty::Float64, Z_matrix::Matrix{Float64})
        return new(environment, chromosome, fitness, objective_function, penalty, Z_matrix)
    end

    function Individual(environment::Environment, chromosome::Vector{Int64}, parent::Individual, new_ones::Vector{Int64}, new_zeros::Vector{Int64})
        execution_statistics["efficient_of_calc_calls"][execution_statistics["current_generation"]] += 1
        if length(new_ones) != length(new_zeros)
            throw(error("Attempt to create Individual with length(ones_changed) != length(zeros_changed). This might indicate that the current solution (chromosome) is unfeasible."))
        end

        if !isdefined(parent, :Z_matrix) || isnothing(parent.Z_matrix)
            try
                execution_statistics["parent_zcalc_from_children"][execution_statistics["current_generation"]] += 1
                parent.Z_matrix = inv(environment.A' * spdiagm(parent.chromosome) * environment.A)

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