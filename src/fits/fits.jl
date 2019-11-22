include("../types.jl")

function batch_fits(treatment_df; num_generations::Int64 = 20000, log_frequency::Int64 = 100, replicates_per_treatment::Int64=50, metadata_file="fits_metadata.csv", data_file="fits.csv", base_random_seed::Int64 = 5, dispersal_kernel_type = @ibd_diskern)

    rseedgenerator = MersenneTwister(base_random_seed)

    CSV.write(metadata_file, treatment_df)

    for (iter, row) in enumerate(eachrow(treatment_df))
        k = row.k
        ibd_str = row.ibd_str
        m = row.m
        ipc = row.ipc
        mu = row.mutation_rate
        n_pops = row.n_pops

        treatment_ct = row.treatment

        @printf("Treatment %d: \t (k=%d, ibd_str=%f, m=%f, ipc=%d, npops=%d)\n\t\t", treatment_ct, k, ibd_str, m, ipc, n_pops)

        @showprogress for rep = 1:replicates_per_treatment
            mp::metapop = init_random_metapop(num_populations=n_pops, diskern_type=dispersal_kernel_type, ibd_decay=ibd_str)
            set_mp_total_k(mp, k)

            rs = rand(rseedgenerator, DiscreteUniform(1, 10^10))

            run_fits(
                    mp,
                    n_gen=num_generations,
                    ipc = ipc,
                    migration_rate = m,
                    mutation_rate = mu,
                    log_freq = log_frequency,
                    rseed = rs,
                    treatment = treatment_ct,
                    replicate = rep,
                    k = k,
                    fits_file = data_file,
                    ibd_str = ibd_str
            )
        end
        println("")
    end
end


function run_fits(mp::metapop; n_gen::Int64=1000, ipc::Int64=5, mutation_rate::Float64 = 10^5, migration_rate::Float64=0.01, log_freq::Int64=20, treatment::Int64=1, replicate::Int64=1, rseed::Int64=1, k::Int64=1000, ibd_str=0, fits_file::String="fits.csv")

    fits_instance::fits = fits(mp, mutation_rate=mutation_rate, migration_rate=migration_rate, n_alleles=ipc, log_freq=log_freq, n_gen=n_gen, rseed=rseed)
    init_fits_uniform_ic(fits_instance)

    n_pops = length(mp.populations)

    df = DataFrame(treatment = [], replicate = [], gen = [], jostd = [], gst = [])

    all_pops = collect(1:n_pops)
    for g = 0:n_gen
        if g % log_freq == 0
            state::Array{Float64} = fits_instance.ct_map
            jostd::Float64 = calc_jost_d(state)
            gst::Float64 = calc_gst(state)
            update_df(df, treatment, replicate, g, jostd, gst)
        end
        run_gen(fits_instance)
    end


    # append if not empty
    CSV.write(fits_file, df, append=isfile(fits_file))
end

function init_fits_uniform_ic(instance::fits)
    n_pops::Int64 = instance.n_pops
    eff_pops::Array{Int64} = [round(2*x.k) for x in instance.mp.populations]
    n_alleles::Int64 = instance.n_alleles
    p = fill(1.0/n_alleles, n_alleles)

    for i = 1:n_pops
        x = rand(instance.rng, Multinomial(eff_pops[i], p))
        for j = 1:n_alleles
            instance.ct_map[i,j] = x[j]
        end
    end
end


function get_new_pj(x_j::Array{Float64}, n_al::Int64, eff_pop_size::Int64)
    if (sum(x_j) > 0 && sum(x_j) != NaN)
        normalize!(x_j, 1)
        xj_new = rand(Multinomial(eff_pop_size, x_j))
        return(xj_new)
    end
    return(fill(0, length(x_j)))
end

function draw_from_diskern_row_old(diskernel_row::Array{Float64})
    u = rand(Uniform())
    max_ind = length(diskernel_row)
    ind::Int64 = 1
    s::Float64 = 0.0
    while (s < u && ind < max_ind)
        s += diskernel_row[ind]
        if s > u
            return ind
        end
        ind += 1
    end

    if (diskernel_row[ind] != 0)
        return ind
    else
        return ind-1
    end
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
        pop_to = draw_from_diskern_row_old(diskernel_row)
        @assert pop_to != pop_from
        for h =1:n_haplo
            eff_pop_size = sum(pop_from_cts)
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
end

function extinct_poly_location(state::Array{Float64,2})
    (x,y) = size(state)
    for j = 1:y
        ext::Bool = true
        for i = 1:x
            if state[i,j] > 0
                ext = false
            end
        end
        if ext
            return j
        end
    end

    return 0
end

function random_mutation(state::Array{Float64, 2}, poly_location)
    (n_pops,n_poly) = size(state)

    rand_pop::Int64 = rand(DiscreteUniform(1,n_pops))

    while (sum(state[rand_pop,:]) < 1)
        rand_pop = rand(DiscreteUniform(1,n_pops))
    end

    rand_i::Int64 = rand(DiscreteUniform(1,n_poly))

    while state[rand_pop, rand_i] < 1
        rand_i = rand(DiscreteUniform(1,n_poly))
    end

    state[rand_pop, rand_i] -= 1
    state[rand_pop, poly_location] += 1
end

function run_gen(instance::fits)
    n_p::Int64 = instance.n_pops
    n_al::Int64 = instance.n_alleles
    migration_rate::Float64 = instance.migration_rate
    mu::Float64 = instance.mutation_rate
    # consider 1 possible mutation per generation
    # only add it if there is an polymorphism location that is empty everywhere

    ext_poly_location = extinct_poly_location(instance.ct_map)
    if ext_poly_location > 0
        if rand(Uniform()) < mu
            random_mutation(instance.ct_map, ext_poly_location)
        end
    end
    #  mutate an individual in a random population


    post_drift::Array{Float64} = zeros(n_p, n_al)
    for p = 1:n_p
        eff_pop_size::Int64 = sum(instance.ct_map[p,:])
        new_pj = get_new_pj(instance.ct_map[p, :], n_al, eff_pop_size)
        post_drift[p,:] = new_pj
    end

    instance.ct_map = post_drift

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

        migration(p, post_drift[p,:], n_indivs_leaving, instance)
    end
end

function update_fits_metadata(metadata::DataFrame, id_ct::Int64, m::Float64, k::Int64, init_poly_ct::Int64, n_pops::Int64)
    push!(metadata.id, id_ct)
    push!(metadata.m, m)
    push!(metadata.k, k)
    push!(metadata.init_poly_ct, init_poly_ct)
    push!(metadata.n_pops, n_pops)
end

function update_df(df::DataFrame, treatment::Int64, replicate::Int64, gen::Int64, jost_d::Float64, gst::Float64)
    push!(df.treatment,treatment)
    push!(df.replicate, replicate)
    push!(df.gen, gen)
    push!(df.jostd, jost_d)
    push!(df.gst, gst)
end


function batch_fits_multicore(treatment_df, n_cores::Int64; num_generations::Int64 = 20000, log_frequency::Int64 = 100, replicates_per_treatment::Int64=50, metadata_file="fits_metadata.csv", data_file="fits.csv", base_random_seed::Int64 = 5, dispersal_kernel_type = @ibd_diskern)

    # split treatment df up according to how many
    treatments_per_core::Int64 = ceil(size(treatment_df)[1] / n_cores)
    @printf("treats: %d, cores: %d, treatments_per_core: %d\n", size(treatment_df)[1], n_cores, treatments_per_core)

    treatment_dfs::Array{DataFrame} = []

    filenames::Array{String} = []
    base_str::String = split(data_file, ".")[1]
    print(base_str)

    for core = 1:n_cores
        filename = string(base_str, "_core", core)
        push!(filenames, filename)


        lo = (core-1)*treatments_per_core + 1
        hi = (core)*treatments_per_core

        if (hi < size(treatment_df)[1])
            hi = size(treatment_df)[1]
        end

        push!(treatment_dfs, treatment_df[lo:hi, :])
    end


    # filename for each treatment

    procs = []

    for (core, df) in enumerate(treatment_dfs)
        base_str = (filenames[core])
        this_metadata::String = string(base_str, "_metadata.csv")
        this_data::String = string(base_str, ".csv")

        s = @spawn batch_fits(
                        df,
                        num_generations=num_generations,
                        replicates_per_treatment = replicates_per_treatment,
                        log_frequency = log_frequency,
                        base_random_seed = 5,
                        metadata_file=this_metadata,
                        data_file=this_data,
                        dispersal_kernel_type=dispersal_kernel_type
                    )
        push!(procs, s)
    end

    for proc in procs
        fetch(proc)
        println(proc)
    end
end
