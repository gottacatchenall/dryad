include("../src/include.jl")

function run_bg_selection()
    m_vals::Array{Float64} = [0.001, 0.005, 0.01]
    s_vals::Array{Float64} = [1.0, 0.5, 0.1]
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
    ibm_metadata::DataFrame = DataFrame(id=[], m=[], s=[], k=[], n_ef=[], n_chromo=[], genome_length=[])
    fits_metadata::DataFrame = DataFrame(id=[], m=[], k=[],init_poly_ct=[])

    lf::Int64 = 20

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
                                    g::genome = init_random_genome(n_ef=n_ef, n_chromo=n_chromo, genome_length=genome_length, init_poly_ct=ipc)


                                    # Run IBM
                                    ibm_instance::ibm = init_ibm(mp, g)
                                    update_ibm_metadata(ibm_metadata, id_ct, m, s, k, n_ef, n_chromo, genome_length)
                                    run_n_generations(ibm_instance, n_gen, migration_rate=m, s=s, genome_file=ibm_genome_file, pop_file=ibm_pop_file, id=id_ct, log_freq=lf)

                                    # Run FITS
                                    n_al::Int64 = convert(Int64, ipc)
                                    fits_instance::fits = fits(mp, n_alleles=convert(Int64,ipc), migration_rate=m, log_freq=lf, n_gen=n_gen)
                                    init_fits_uniform_ic(fits_instance)

                                    update_fits_metadata(fits_metadata, id_ct, m, k, ipc)
                                    df = run_fits(fits_instance, id_ct)
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
    CSV.write("bg_ibm_metadata.csv", ibm_metadata)
    CSV.write("bg_fits_metadata.csv", fits_metadata)
end

run_bg_selection()
