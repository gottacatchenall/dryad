include("./src/include.jl")

function run_test(; n_gen::Int64=500)
    mp::metapop = init_random_metapop()
    g::genome = init_random_genome()
    instance::ibm = init_ibm(mp, g)

    @time run_n_generations(instance, n_gen)
    #Profile.print()
end

function run_batch_ibm()
    m_vals::Array{Float64} = [0.005]
    s_vals::Array{Float64} = [3.0, 5.0, 8.0]
    k_vals::Array{Int64} = [500]
    init_poly_ct::Array{Float64} = [3.0, 8.0, 15.0]
    n_ef_vals::Array{Int64} = [1]
    n_chromo_vals::Array{Int64} = [5]
    genome_length_vals::Array{Float64} = [100.0]
    n_rep = 1
    n_gen = 1000

    id_ct::Int64 = 1

    ibm_pop_file::String = "ibm_pops.csv"
    ibm_genome_file::String = "ibm_genomes.csv"
    fits_file::String = "fits.csv"

    df = DataFrame(id=[],gen=[],jostd=[],gst=[])
    CSV.write(fits_file, df)


    # need to track metadata
    ibm_metadata::DataFrame = DataFrame(id=[], m=[], s=[], k=[], n_ef=[], n_chromo=[], genome_length=[], init_poly_ct=[])
    fits_metadata::DataFrame = DataFrame(id=[], m=[], k=[],init_poly_ct=[])

    lf::Int64 = 20

    # init with appropriate ef values?

    base_random_seed = 5
    rseedgenerator = MersenneTwister(base_random_seed)

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

                                    rs = rand(rseedgenerator, DiscreteUniform(1, 10^10))

                                    ## ================================
                                    ## Run IBM
                                    ## ================================
                                    println("ID: ", id_ct, "\t m=",m," k=",k, " s=", s, " ipc=", ipc)
                                    print("\tIBM:  ", )
                                    @time run_ibm(mp, ibm_metadata, n_gen=n_gen, migration_rate=m, selection_strength = s, log_freq=lf, genome_file=ibm_genome_file, pop_file=ibm_pop_file, id=id_ct, init_poly_ct=ipc, k=k, rseed=rs)


                                    ## ================================
                                    ## Run FITS
                                    ## ================================
                                    print("\tFITS:  ")
                                    # function run_fits(mp::metapop; n_gen::Int64=1000, ipc::Int64=5, migration_rate::Float64=0.01, log_freq::Int64=20)
                                    @time run_fits(mp, fits_metadata, n_gen=n_gen, ipc=convert(Int64,ipc), migration_rate=m, log_freq=lf, rseed=rs, id=id_ct, k=k, fits_file=fits_file)


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

function fitstest()
    n_als = (3, 8, 15)
    n_pops = (20)
    m = (0.001,0.005, 0.01)
    eff_pop_size_list = (200, 800, 2000)
    ng = 20000
    n_rep = 30

    fits_metadata::DataFrame = DataFrame(id=[], m=[], k=[],init_poly_ct=[])


    df = DataFrame()
    df.id = []
    df.gen = []
    df.jostd = []
    df.gst = []
    CSV.write("output.csv", df)

    lf = 1000
    idct::Int64 = 0
    mp::metapop = init_random_metapop()
    for base_mig in m
        for n_al in n_als
            for n_pop in n_pops
                for eff_pop in eff_pop_size_list
                    set_mp_total_k(mp, eff_pop)
                    for rep = 1:n_rep

                        fits_instance::fits = fits(mp, migration_rate=base_mig, log_freq=lf, n_gen=ng, n_alleles=n_al)
                        init_fits_uniform_ic(fits_instance)
                        #df = run_dke(n_gen, n_pop, n_al, base_mig, eff_pop, 10)
                        df = run_fits(fits_instance, idct)

                        update_fits_metadata(fits_metadata, idct, base_mig, eff_pop, convert(Float64,n_al))


                        idct += 1
                        CSV.write("output.csv", df, append=true)
                    end
                end
            end
        end
    end
    CSV.write("metadata.csv", fits_metadata)
end

#@time fitstest()
@time run_batch_ibm()
#run_test()

#run_fits_test()
