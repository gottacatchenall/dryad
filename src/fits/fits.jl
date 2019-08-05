include("../types.jl")

function get_new_pj(x_j::Array{Float64}, n_al::Int64, eff_pop_size::Int64)
    if (sum(x_j) > 0 )
        normalize!(x_j, 1)
        xj_new = rand(Multinomial(eff_pop_size, x_j))
        return(xj_new)
    else
        return x_j
    end
end

function init_fits_uniform_ic(instance::fits)
    n_pops::Int64 = instance.n_pops
    eff_pops::Array{Int64} = [round(2*x.k) for x in instance.mp.populations]
    n_alleles::Int64 = instance.n_alleles

    for p = 1:n_pops
        base_ct::Int64 = round(eff_pops[p]/n_alleles)
        for i = 1:n_alleles
            instance.ct_map[p,i] = base_ct
        end
    end
    # find a way to do divison, or something
end


function run_fits(instance::fits, log_freq::Int64, n_gen::Int64, id::Int64)
    df = DataFrame(id=[],gen=[],gst=[],jostd=[])
    n_pops::Int64 = length(instance.mp.populations)

    for g = 0:n_gen
        if g % log_freq == 0
            state::Array{Float64} = instance.ct_map
            jostd::Float64 = calc_jost_d(state)
            gst::Float64 = calc_gst(state)
            push!(df.jostd, jostd)
            push!(df.gen, g)
            push!(df.gst, gst)
            push!(df.id, id)
        end
        run_gen(instance)
    end

    return df
end

function get_allele_list(pop_alleles::Array{Float64}, eff_pop_size::Int64)
    allele_list::Array{Int64} = zeros(eff_pop_size)

    n_alleles::Int64 = length(pop_alleles)

    ct::Int64 = 1
    for i = 1:n_alleles
        new = pop_alleles[i]
        for j = 1:new
            allele_list[ct] = i
            ct += 1
        end
    end
    return(allele_list)
end

function migration(pop_from::Int64, pop_from_cts::Array{Float64}, n_indivs_leaving::Int64, instance::fits)
    diskernel_row::Array{Float64} = instance.mp.diskern.D[pop_from,:]
    n_pops::Int64 = instance.n_pops
    n_alleles::Int64 = instance.n_alleles
    n_haplo::Int64 = 2
    eff_pop_size::Int64 = sum(pop_from_cts)

    if (eff_pop_size/2 < n_indivs_leaving)
        n_indivs_leaving = eff_pop_size/2
    end


    for ind = 1:n_indivs_leaving
        #draw_from_diskern_row(pop::Int64, row::Array{Float64,1}, n_pops::Int64)

        pop_to = draw_from_diskern_row(pop_from, diskernel_row, n_pops)
        @assert pop_to != pop_from
        for h =1:n_haplo
            eff_pop_size = sum(pop_from_cts)

            # ya this isn't the best data structure to use but
            allele_list::Array{Int64} = get_allele_list(pop_from_cts, eff_pop_size)
            random_index = rand(DiscreteUniform(1, eff_pop_size))

            allele_num = allele_list[random_index]

            pop_from_cts[allele_num] -= 1
            instance.ct_map[pop_from, allele_num] -= 1
            instance.ct_map[pop_to, allele_num] += 1

            @assert instance.ct_map[pop_from, allele_num] >= 0
            @assert instance.ct_map[pop_to, allele_num] >= 0
        end
    end

    return(instance)
end

function run_gen(instance::fits)
    n_p::Int64 = instance.n_pops
    n_al::Int64 = instance.n_alleles
    migration_rate::Float64 = instance.migration_rate

    # drift
    post_drift::Array{Float64} = zeros(n_p, n_al)

    # store this in a new array
    for p = 1:n_p
        eff_pop_size::Int64 = sum(instance.ct_map[p,:])
        new_pj = get_new_pj(instance.ct_map[p, :], n_al, eff_pop_size)
        post_drift[p,:] = new_pj
    end

    #post_drift_now = deepcopy(post_drift)
    instance.ct_map = post_drift
    # migration

    # make this a binomial you spoon
    for p = 1:n_p
        post_drift_this_pop = (post_drift[p,:])
        eff_pop_size::Int64 = sum(post_drift_this_pop)
        exp_n_alleles_leaving::Float64 = (eff_pop_size/2.0) * migration_rate
        rem::Float64 = exp_n_alleles_leaving - floor(exp_n_alleles_leaving)
        extra_mig::Int64 = 0
        if (rand(Uniform()) < rem)
            extra_mig = 1
        end
        n_indivs_leaving::Int64 = (floor(exp_n_alleles_leaving) + extra_mig)

        # pick random in 1...eff_pop, and that in the ct is the allele
        # update the real df, pass it the post drift matrix to compute where migrants go,

        migration(p, post_drift[p,:], n_indivs_leaving, instance)
    end
end
