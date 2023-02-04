mutable struct Individual
    ga::GeneticAlgorithm
    chromosome::Vector{Int64}
    fitness::Float64
    objective_function::Float64
    penalty::Float64
    

    function Individual(ga::GeneticAlgorithm, chromosome::Vector{Int64})
        fitness, objective_function, penalty = calculate_fitness(chromosome, ga.A)
        new(ga, chromosome, fitness, objective_function, penalty)
    end

end

function ldet(A)
    (value,sign_det) = logabsdet(A);
    if (sign_det > 0)
        return value
    else
        return -Inf
    end
end

function calculate_fitness(chromosome, A)
    num_ones = sum(chromosome)
    objective_function = ldet(A'*diagm(vec(chromosome))*A);
    penalty = - 100 * abs(num_ones - 5)
    fitness = objective_function + penalty
    return fitness, objective_function, penalty
end

function mutate(ga::GeneticAlgorithm, individual::Individual)
    
end

function breed(parent1::Individual, parent2::Individual)

end