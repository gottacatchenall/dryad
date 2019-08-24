include("../src/include.jl")

function run_bg_selection()
    m_vals::Array{Float64} = [0.001, 0.005, 0.01]
    s_vals::Array{Float64} = [3.0,5.0,8.0]
    k_vals::Array{Int64} = [500,2000,4000]
    init_poly_ct::Array{Float64} = [3.0, 8.0, 15.0]
    n_ef_vals::Array{Int64} = [1]
    n_chromo_vals::Array{Int64} = [5]
    genome_length_vals::Array{Float64} = [100.0]
    n_rep = 50
    n_gen = 1000

    id_ct::Int64 = 1

    ibm_pop_file::String = "bg_ibm_pops.csv"
    ibm_genome_file::String = "bg_ibm_genomes.csv"
    fits_file::String = "bg_fits.csv"

    df = DataFrame(id=[],gen=[],jostd=[],gst=[])
    CSV.write(fits_file, df)


    # need to track metadata
    ibm_metadata::DataFrame = DataFrame(id=[], m=[], s=[], k=[], n_ef=[], n_chromo=[], genome_length=[], init_poly_ct=[])
    fits_metadata::DataFrame = DataFrame(id=[], m=[], k=[],init_poly_ct=[])

    lf::Int64 = 20

    base_random_seed = 5
    rseedgenerator = MersenneTwister(base_random_seed)

    # init with appropriate ef values?

    mp::metapop = init_random_metapop(selection_type=@background_selection)
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
                                    println("ID: ", id_ct, "\t\t\t m=",m," k=",k, " s=", s, " ipc=", ipc, " rs=", rs)
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
    CSV.write("bg_ibm_metadata.csv", ibm_metadata)
    CSV.write("bg_fits_metadata.csv", fits_metadata)
end

run_bg_selection()
