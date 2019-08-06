include("./src/include.jl")

function run_test(; n_gen::Int64=500)
    mp::metapop = init_random_metapop()
    g::genome = init_random_genome()
    instance::ibm = init_ibm(mp, g)

    @time run_n_generations(instance, n_gen)
    #Profile.print()
end

function update_ibm_metadata(metadata::DataFrame,id_ct::Int64, m::Float64, s::Float64, k::Int64, n_ef::Int64, n_chromo::Int64, genome_length::Float64)
    push!(metadata.id, id_ct)
    push!(metadata.m, m)
    push!(metadata.s, s)
    push!(metadata.k, k)
    push!(metadata.n_ef, n_ef)
    push!(metadata.n_chromo, n_chromo)
    push!(metadata.genome_length, genome_length)
end

function update_fits_metadata(metadata::DataFrame, id_ct::Int64, m::Float64, k::Int64)
    push!(metadata.id, id_ct)
    push!(metadata.m, m)
    push!(metadata.k, k)
end

function run_batch_ibm()
    m_vals::Array{Float64} = [0.001, 0.005, 0.01]
    s_vals::Array{Float64} = [1.0]
    k_vals::Array{Int64} = [500, 2000, 4000]
    init_poly_ct::Array{Float64} = [3.0, 8.0, 15.0]
    n_ef_vals::Array{Int64} = [1]
    n_chromo_vals::Array{Int64} = [5]
    genome_length_vals::Array{Float64} = [100.0]
    n_rep = 1
    n_gen = 500

    id_ct::Int64 = 1

    ibm_pop_file::String = "ibm_pops.csv"
    ibm_genome_file::String = "ibm_genomes.csv"
    fits_file::String = "fits.csv"

    df = DataFrame(id=[],gen=[],gst=[],jostd=[])
    CSV.write(fits_file, df)


    # need to track metadata
    ibm_metadata::DataFrame = DataFrame(id=[], m=[], s=[], k=[], n_ef=[], n_chromo=[], genome_length=[])
    fits_metadata::DataFrame = DataFrame(id=[], m=[], k=[])

    lf::Int64 = 20

    mp::metapop = init_random_metapop()
    for ipc in init_poly_ct
        for m in m_vals
            for s in s_vals
                for k in k_vals
                    for n_ef in n_ef_vals
                        for n_chromo in n_chromo_vals
                            for genome_length in genome_length_vals
                                for r = 1:n_rep
                                    set_mp_total_k(mp, k)
                                    g::genome = init_random_genome(n_ef=n_ef, n_chromo=n_chromo, genome_length=genome_length, init_poly_ct=ipc)


                                    # Run IBM
                                    ibm_instance::ibm = init_ibm(mp, g)
                                    update_ibm_metadata(ibm_metadata, id_ct, m, s, k, n_ef, n_chromo, genome_length)
                                    run_n_generations(ibm_instance, n_gen, migration_rate=m, s=s, genome_file=ibm_genome_file, pop_file=ibm_pop_file, id=id_ct, log_freq=lf)

                                    # Run FITS
                                    n_al::Int64 = convert(Int64, ipc)
                                    fits_instance::fits = fits(mp, n_al, m)
                                    init_fits_uniform_ic(fits_instance)

                                    update_fits_metadata(fits_metadata, id_ct, m, k)
                                    df = run_fits(fits_instance, lf, n_gen, id_ct)
                                    CSV.write(fits_file, df, append=true)

                                    id_ct += 1
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    CSV.write("ibm_metadata.csv", ibm_metadata)
    CSV.write("fits_metadata.csv", fits_metadata)
end

function set_mp_total_k(mp::metapop, k::Int64)
    n_pops::Int64 = length(mp.populations)
    base::Float64 = (k/n_pops)
    for p = 1:n_pops
        mp.populations[p].k = base
    end
    mp.n_indivs = k
end

function run_fits_test()
    mp::metapop = init_random_metapop()

    n_pops = length(mp.populations)
    n_gen::Int64 = 500
    n_alleles::Int64 = 5
    m::Float64 = 0.01
    eff_pop::Float64 = 20.0
    log_freq::Int64 = 20

#    fits(n_pops::Int64, n_alleles::Int64, migration_rate::Float64, eff_pop_size::Int64, mp::metapop) = new(n_pops, n_alleles, migration_rate, eff_pop_size, zeros(Float64, n_pops, n_alleles), mp)

    instance::fits = fits(mp, n_alleles, m)

    # gotta initialize it with some type of state, or somthing
    init_fits_uniform_ic(instance)

    df = run_fits(instance, log_freq, 300)
    CSV.write("fitstest.csv", df)
end

run_batch_ibm()
#run_test()

#run_fits_test()
