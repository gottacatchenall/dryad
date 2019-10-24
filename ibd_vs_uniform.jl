include("./src/include.jl")

function fits()
    n_als = (10)
    n_pops = (20)
    mig_rates = collect(0.0001:0.0001:0.01)
    k_vals = (2000)
    n_gen = 20000
    n_rep = 50
    lf = 500

    base_random_seed = 5
    rseedgenerator = MersenneTwister(base_random_seed)



    ## ===============================================================
    #
    #           IBD dispersal kernel
    #
    ## ===============================================================


    fits_metadata::DataFrame = DataFrame(id=[], m=[], k=[],init_poly_ct=[],n_pops=[])
    fits_file::String = "ibd_fits.csv"
    df = DataFrame(id=[],gen=[],jostd=[],gst=[])
    CSV.write(fits_file, df)
    idct::Int64 = 0
    for k in k_vals
        for m in mig_rates
            for ipc in n_als
                for n_pop in n_pops
                    mp::metapop = init_random_metapop(num_populations=n_pop, diskern_type=@ibd_diskern)
                    set_mp_total_k(mp, k)
                    for rep = 1:n_rep
                        println("FITS \t\tID: ", idct, "\t m=",m," k=",k, " ipc=", ipc)
                        print("\t\t")
                        rs = rand(rseedgenerator, DiscreteUniform(1, 10^10))
                        @time run_fits(mp, fits_metadata, n_gen=n_gen, ipc=convert(Int64,ipc), migration_rate=m, log_freq=lf, rseed=rs, id=idct, k=k, fits_file=fits_file)
                        idct = idct + 1
                    end
                end
            end
        end
    end
    CSV.write("ibd_metadata.csv", fits_metadata)


    ## ===============================================================
    #
    #           uniform dispersal kernel
    #
    ## ===============================================================

    fits_metadata = DataFrame(id=[], m=[], k=[],init_poly_ct=[],n_pops=[])
    fits_file = "uniform_fits.csv"
    df = DataFrame(id=[],gen=[],jostd=[],gst=[])
    CSV.write(fits_file, df)
    for k in k_vals
        for m in mig_rates
            for ipc in n_als
                for n_pop in n_pops
                    mp::metapop = init_random_metapop(num_populations=n_pop, diskern_type=@uniform_diskern)
                    set_mp_total_k(mp, k)
                    for rep = 1:n_rep
                        println("FITS \t\tID: ", idct, "\t m=",m," k=",k, " ipc=", ipc)
                        print("\t\t")
                        rs = rand(rseedgenerator, DiscreteUniform(1, 10^10))
                        @time run_fits(mp, fits_metadata, n_gen=n_gen, ipc=convert(Int64,ipc), migration_rate=m, log_freq=lf, rseed=rs, id=idct, k=k, fits_file=fits_file)
                        idct = idct + 1
                    end
                end
            end
        end
    end
    CSV.write("uniform_metadata.csv", fits_metadata)
end

fits()
