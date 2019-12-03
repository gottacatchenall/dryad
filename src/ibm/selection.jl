
function selection(instance::ibm; b::Float64=3.0, s::Float64=3.0)
    pop_cts::Array{Float64, 1} = get_cts(instance)
    n_indivs::Int64 = length(instance.population_map)
    n_pops::Int64 = length(pop_cts)
    pop::Int64 = 0

    for i = 1:n_indivs
        pop = instance.population_map[i]
        if (pop != 0)
            k::Float64 = instance.mp.populations[pop].k
            n::Float64 = pop_cts[pop]
            efs::Array{Float64} = instance.mp.populations[pop].efs
            w::Float64 = calc_fitness(i, efs, instance, s)
            instance.fitness_map[i] = w

            surv::Bool = beverton_holt(w, n, k, b)
            if (!surv)
                instance.population_map[i] = 0
            end
        end
    end
end

function beverton_holt(w::Float64, n::Float64, k::Float64, b::Float64)
    prob::Float64 = 1.0/(1.0 + (b - 1)*(n/k))
    u = rand()
    if (u < prob)
        return true
    else
        return false
    end
end

function calc_fitness(i::Int64, efs::Array{Float64}, instance::ibm, s::Float64)::Float64
    chromos::Array{chromosome} = instance.g.chromosomes
    n_chromo::Int64 = length(chromos)
    genotypes::Array{Float64, 3} = instance.genotypes

    # normalize total fitness by n fitness loci


    w::Float64 = 1.0
    l::Int64 = 1
    for c = 1:n_chromo
        this_chr_n_loci::Int64 = chromos[c].n_loci
        for c_l = 1:this_chr_n_loci
            corresponding_ef::Int64 = chromos[c].ef_map[c_l]
            #println("corresponding_ef: ", corresponding_ef)
            if (corresponding_ef != 0)
                ef_val = efs[corresponding_ef]

                for h = 1:2
                    x_i::Float64 = genotypes[i, l, h]
                    w_i = calc_fitness_component(x_i, ef_val, s)
                    w = w * w_i
                end
            end
            l += 1
        end
    end
    #println("indiv, fitness: ", i, " ", w)

    return w
end

function calc_fitness_component(x_i::Float64, ef_val::Float64, s::Float64)
    diff_sq::Float64 = (abs(x_i - ef_val))^2

#    println("x_i, ef, diff_sq: ", x_i, " ", ef_val," ", diff_sq)

    s_inv::Float64 = 1.0 / (2*(s^2))

    gauss::Float64 = exp((-1*(diff_sq*s_inv)))
    sel_str::Float64 = 1.0
    return ((sel_str*gauss))
end

function get_cts(instance::ibm)::Array{Float64,1}
    n_pops::Int64 = length(instance.mp.populations)
    prop_full::Array{Float64,1} = zeros(n_pops)

    max_i::Int64 = length(instance.population_map)
    pop::Int64 = 0

    for i = 1:max_i
        pop = instance.population_map[i]
        if (pop > 0)
            prop_full[pop] += 1
        end
    end
    return prop_full
end
