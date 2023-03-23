
function reset_commoner_population(environment::Environment, elite::Vector{Individual})
    # println("Resetting commoners.")
    commoners = initialize_population(environment) # could temporarily change population size to avoid generating more individuals than needed
    needed_commoners = environment.population_size - length(elite)
    # println("needed_commoners ", needed_commoners)
    return [elite; commoners[1:needed_commoners]]
end

function adapt(environment::Environment, elite::Vector{Individual})
    population = Vector{Individual}()

    if environment.adaptation_method == "reset"
        population = reset_commoner_population(environment, elite)
    else
        throw(error("Unknow adaptation method ", environment.adaptation_method))
    end

    return population
end