
function dispersal(instance::ibm; migration_rate::Float64=0.01)
    current_state::Array{Int64, 1} = copy(instance.population_map)
    n_pops::Int64 = length(instance.mp.populations)
    max_n_indivs = length(current_state)

    indivs_by_pop::Array{Array{Int64,1}} = split_by_pop(current_state, n_pops)

    for p = 1:n_pops
        n_indivs::Int64 = length(indivs_by_pop[p])
        n_mig::Int64 = get_number_migrants_this_pop(n_indivs, migration_rate)

        for ind = 1:n_mig
            new_pop::Int64 = draw_from_diskern_row(p, instance.dk.D[p,:], n_pops)
            mig_index::Int64 = get_migrant(indivs_by_pop[p])
            instance.population_map[mig_index] = new_pop
        end
    end
end

function split_by_pop(current_state::Array{Int64,1}, n_pops::Int64)::Array{Array{Int64 ,1}}
    indivs_by_pop::Array{Array{Int64,1}} = fill([], n_pops)
    n::Int64 = length(current_state)
    this_indiv_pop::Int64 = 0
    for i = 1:n
        this_indiv_pop = current_state[i]
        if (this_indiv_pop > 0)
            push!(indivs_by_pop[this_indiv_pop], i)
        end
    end

    return indivs_by_pop
end

function get_migrant(current_pop::Array{Int64, 1})
    # only pull from indecies where current_pop is non zero
    mig_index_this_pop::Int64 = 0
    mig_index_global::Int64 = 0
    while (mig_index_global == 0)
        mig_index_this_pop = rand(1:length(current_pop))
        mig_index_global = current_pop[mig_index_this_pop]
    end
    current_pop[mig_index_this_pop] = 0
    return mig_index_global
end

function get_number_migrants_this_pop(n_indivs::Int64, migration_rate::Float64)
    return rand(Binomial(n_indivs, migration_rate))
end


function draw_from_diskern_row(pop::Int64, row::Array{Float64,1}, n_pops::Int64)
    u::Float64 = rand(Uniform())

    s::Float64 = 0.0
    for p2 = 1:n_pops
        if p2 != pop
            s += row[p2]
            if (s > u)
                return p2
            end
        end
    end
    if (pop != n_pops)
        return (n_pops)
    else
        return(1)
    end

end
